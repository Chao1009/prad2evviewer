//=============================================================================
// HyCalCluster.cpp — Island clustering for HyCal
//
// Ported from PRadIslandCluster / PRadHyCalReconstructor (PRadAnalyzer).
// Uses pre-computed neighbor lists from HyCalSystem for fast adjacency checks.
//
// Chao Peng (original PRadAnalyzer), adapted for prad2decoder.
//=============================================================================

#include "HyCalCluster.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace fdec
{

static constexpr int ISLAND_GROUP_RESERVE = 50;
static constexpr int POS_RECON_HITS       = 15;

namespace
{

// Log-weighted centroid around a centre module.  Offsets are in module sizes;
// the centre itself, at (0, 0), seeds the sums with its own weight.
struct LogCentroid {
    float wx = 0.f, wy = 0.f, wtot;
    int   n;

    explicit LogCentroid(float center_weight)
        : wtot(center_weight), n(center_weight > 0.f ? 1 : 0) {}

    void add(float dx, float dy, float w)
    {
        if (w > 0.f) {
            wx += dx * w;
            wy += dy * w;
            wtot += w;
            ++n;
        }
    }

    // position in mm; the module centre when nothing carries weight
    void position(const Module &center, float &x, float &y) const
    {
        if (wtot > 0.f) {
            x = center.x + (wx / wtot) * center.size_x;
            y = center.y + (wy / wtot) * center.size_y;
        } else {
            x = center.x;
            y = center.y;
        }
    }
};

} // namespace

float shower_depth(int center_id, float energy_mev)
{
    if (energy_mev <= 0.f) return 0.f;
    if (center_id >= PWO_ID0)            // PbWO4 (W-modules)
        return 8.6f  * (std::log(energy_mev / 1.1f)  - 0.5f);
    return            26.7f * (std::log(energy_mev / 2.84f) - 0.5f);  // PbGlass
}

HyCalCluster::HyCalCluster(const HyCalSystem &sys)
    : sys_(sys)
    , profile_(std::make_shared<SimpleProfile>())
{
}

HyCalCluster::~HyCalCluster() = default;

void HyCalCluster::SetProfile(std::shared_ptr<const IClusterProfile> prof)
{
    profile_ = prof ? std::move(prof) : std::make_shared<SimpleProfile>();
}

void HyCalCluster::SetProfile(IClusterProfile *prof)
{
    SetProfile(std::shared_ptr<const IClusterProfile>(prof));
}

// --- per-event interface ----------------------------------------------------

void HyCalCluster::Clear()
{
    hits_.clear();
    groups_.clear();
    clusters_.clear();
}

void HyCalCluster::AddHit(int module_index, float energy, float time)
{
    if (module_index < 0 || module_index >= sys_.module_count()) return;
    if (test_bit(sys_.module(module_index).flag, kDeadModule)) return;
    if (energy > config_.min_module_energy)
        hits_.push_back({module_index, energy, time});
}

void HyCalCluster::FormClusters()
{
    clusters_.clear();
    groups_.clear();

    group_hits();

    // find maxima and split each group
    for (auto &group : groups_)
        split_cluster(group);

    if (config_.leakage_correction) {
        for (auto &cl : clusters_)
            apply_leakage_correction(cl);
    }
}

void HyCalCluster::ReconstructHits(std::vector<ClusterHit> &out) const
{
    out.clear();
    out.reserve(clusters_.size());

    for (auto &cl : clusters_) {
        if (cl.energy < config_.min_cluster_energy) continue;
        if (static_cast<int>(cl.hits.size()) < config_.min_cluster_size) continue;
        out.push_back(reconstruct_pos(cl));
    }
    std::sort(out.begin(), out.end(), [](const ClusterHit &a, const ClusterHit &b) {
        return a.energy > b.energy;
    });
}

void HyCalCluster::ReconstructMatched(std::vector<RecoResult> &out) const
{
    out.clear();
    out.reserve(clusters_.size());

    for (auto &cl : clusters_) {
        if (cl.energy < config_.min_cluster_energy) continue;
        if (static_cast<int>(cl.hits.size()) < config_.min_cluster_size) continue;
        out.push_back({&cl, reconstruct_pos(cl)});
    }
}

// Seed-driven BFS grouping (multi-pulse aware); the algorithm is described at
// ClusterConfig::seed_time_window in HyCalCluster.h.
void HyCalCluster::group_hits()
{
    // Build per-module pulse lists.  Multiple pulses on the same module
    // are common with FADC waveform data; mod_to_hits_[m] holds every
    // hit index whose ModuleHit::index == m.
    mod_to_hits_.assign(sys_.module_count(), {});
    for (int i = 0; i < static_cast<int>(hits_.size()); ++i)
        mod_to_hits_[hits_[i].index].push_back(i);

    hit_group_id_.assign(hits_.size(), -1);
    consumed_.assign(hits_.size(), false);

    // Energy-descending seed order — the global largest pulse seeds first,
    // ensuring that if multiple pulses on a module are time-coincident with
    // the seed, the most energetic shower claims the cluster.
    std::vector<int> order(hits_.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
              [&](int a, int b) { return hits_[a].energy > hits_[b].energy; });

    for (int seed : order) {
        if (consumed_[seed]) continue;
        // sorted desc — once we drop below the seed threshold no later
        // pulse can satisfy it.
        if (hits_[seed].energy < config_.min_center_energy) break;

        int gid = static_cast<int>(groups_.size());
        groups_.emplace_back();
        groups_.back().reserve(ISLAND_GROUP_RESERVE);

        consumed_[seed] = true;
        hit_group_id_[seed] = gid;
        grow_island(seed, gid, groups_.back());
    }
}

void HyCalCluster::grow_island(int seed_idx, int group_id, std::vector<int> &group)
{
    const float seed_time = hits_[seed_idx].time;
    const bool  use_time  = config_.seed_time_window > 0.f;
    const float dt_max    = config_.seed_time_window;

    // BFS using the group vector itself as the queue: indices [qi..size())
    // are still to be expanded.  Seed pushed by caller; we expand the
    // frontier in insertion order.
    group.push_back(seed_idx);
    for (size_t qi = 0; qi < group.size(); ++qi) {
        int hi = group[qi];
        sys_.for_each_neighbor(hits_[hi].index, config_.corner_conn, [&](int ni) {
            // Per neighbour MODULE, pick at most one pulse to add to this
            // group: the LARGEST-energy unconsumed pulse whose time lies
            // within ±seed_time_window of the seed (or any unconsumed
            // pulse when gating is off).  Other pulses on the same module
            // remain in the pool for a different seed at a different timing.
            int best_k = -1;
            for (int k : mod_to_hits_[ni]) {
                if (consumed_[k]) continue;
                if (use_time && std::fabs(hits_[k].time - seed_time) > dt_max)
                    continue;
                if (best_k < 0 || hits_[k].energy > hits_[best_k].energy)
                    best_k = k;
            }
            if (best_k >= 0) {
                consumed_[best_k] = true;
                hit_group_id_[best_k] = group_id;
                group.push_back(best_k);
            }
        });
    }
}

// --- split cluster: find local maxima, distribute hits ----------------------

void HyCalCluster::split_cluster(const std::vector<int> &group)
{
    auto maxima = find_maxima(group);
    if (maxima.empty()) return;

    if (maxima.size() == 1 ||
        group.size() >= SPLIT_MAX_HITS ||
        maxima.size() >= SPLIT_MAX_MAXIMA)
    {
        // single cluster from this group
        auto &seed = hits_[maxima[0]];
        clusters_.emplace_back();
        auto &cl = clusters_.back();
        cl.center = seed;
        cl.flag   = sys_.module(seed.index).flag;

        for (int hi : group)
            cl.add_hit(hits_[hi]);
        cl.energy_square = calculate_energy_square(cl);
    }
    else {
        split_hits(maxima, group);
    }
}

std::vector<int> HyCalCluster::find_maxima(const std::vector<int> &group) const
{
    std::vector<int> local_max;
    local_max.reserve(20);
    if (group.empty()) return local_max;

    const int gid = hit_group_id_[group[0]];

    for (int hi : group) {
        auto &hit = hits_[hi];
        if (hit.energy < config_.min_center_energy)
            continue;

        bool is_max = true;
        // include corners when checking for maxima.  With multi-pulse
        // modules, only the pulse that joined this group counts — others on
        // the same module belong to a different (later) seed.
        sys_.for_each_neighbor(hit.index, true, [&](int ni) {
            if (!is_max) return;
            for (int hj : mod_to_hits_[ni]) {
                if (hit_group_id_[hj] == gid && hits_[hj].energy > hit.energy) {
                    is_max = false;
                    return;
                }
            }
        });

        if (is_max)
            local_max.push_back(hi);
    }

    return local_max;
}

// Distribute shared hits among multiple maxima.
void HyCalCluster::split_hits(const std::vector<int> &maxima,
                               const std::vector<int> &group)
{
    SplitContainer split;  // ~4KB on stack, safe per-call

    int nmax  = static_cast<int>(maxima.size());
    int nhits = static_cast<int>(group.size());

    // initialize fractions from profile
    for (int i = 0; i < nmax; ++i) {
        auto &center = hits_[maxima[i]];
        for (int j = 0; j < nhits; ++j) {
            auto &hit = hits_[group[j]];
            split.frac[j][i] = get_profile_frac(center, hit) * center.energy;
        }
    }

    // iterative refinement
    eval_fraction(maxima, group, split);

    // create clusters from final fractions
    for (int i = 0; i < nmax; ++i) {
        clusters_.emplace_back();
        auto &cl = clusters_.back();
        cl.center = hits_[maxima[i]];
        cl.flag   = sys_.module(cl.center.index).flag;

        for (int j = 0; j < nhits; ++j) {
            if (split.frac[j][i] == 0.f) continue;

            float nf = split.norm_frac(i, j);
            if (nf < config_.least_split) {
                split.total[j] -= split.frac[j][i];
                continue;
            }

            ModuleHit new_hit = hits_[group[j]];
            new_hit.energy *= nf;
            cl.add_hit(new_hit);

            if (new_hit.index == cl.center.index)
                cl.center.energy = new_hit.energy;

            set_bit(cl.flag, kSplit);
        }
        cl.energy_square = calculate_energy_square(cl);
    }
}

float HyCalCluster::calculate_energy_square(const ModuleCluster &cluster) const
{
    const auto &center_mod = sys_.module(cluster.center.index);
    float energy_square = 0.f;

    for (const auto &hit : cluster.hits) {
        if (config_.seed_time_window > 0.f &&
            std::fabs(hit.time - cluster.center.time) > config_.seed_time_window)
            continue;

        double dx, dy;
        sys_.qdist(center_mod, sys_.module(hit.index), dx, dy);
        if (qdist_in_5x5(dx, dy))
            energy_square += hit.energy;
    }
    return energy_square;
}

void HyCalCluster::eval_fraction(const std::vector<int> &maxima,
                                  const std::vector<int> &group,
                                  SplitContainer &split) const
{
    int nmax  = static_cast<int>(maxima.size());
    int nhits = static_cast<int>(group.size());

    struct BaseHit { float x, y, E; };
    BaseHit temp[POS_RECON_HITS];

    int iters = config_.split_iter;
    while (iters-- > 0) {
        split.sum_frac(nhits, nmax);

        for (int i = 0; i < nmax; ++i) {
            auto &center_hit = hits_[maxima[i]];
            const auto &center_mod = sys_.module(center_hit.index);

            // gather 3x3 neighbors for position reconstruction
            float tot_E = center_hit.energy;
            int count = 0;

            for (int j = 0; j < nhits; ++j) {
                auto &hit = hits_[group[j]];
                if (hit.index == center_hit.index || split.frac[j][i] == 0.f)
                    continue;

                const auto &hit_mod = sys_.module(hit.index);
                double dx, dy;
                sys_.qdist(center_mod, hit_mod, dx, dy);

                if (qdist_in_3x3(dx, dy) && count < POS_RECON_HITS) {
                    float frac_E = hit.energy * split.norm_frac(i, j);
                    temp[count] = {static_cast<float>(dx), static_cast<float>(dy), frac_E};
                    tot_E += frac_E;
                    count++;
                }
            }

            // reconstruct position (log-weighted)
            LogCentroid acc(get_weight(center_hit.energy, tot_E));
            for (int k = 0; k < count; ++k)
                acc.add(temp[k].x, temp[k].y, get_weight(temp[k].E, tot_E));

            // the centre is rounded to float before the shift, unlike
            // LogCentroid::position; the split fractions depend on it
            float cx = center_mod.x, cy = center_mod.y;
            if (acc.wtot > 0.f) {
                cx += (acc.wx / acc.wtot) * center_mod.size_x;
                cy += (acc.wy / acc.wtot) * center_mod.size_y;
            }

            // update fractions with new center position
            for (int j = 0; j < nhits; ++j) {
                auto &hit = hits_[group[j]];
                split.frac[j][i] = ProfileFractionAt(cx, cy, tot_E, hit.index) * tot_E;
            }
        }
    }
    split.sum_frac(nhits, nmax);
}

// --- position reconstruction: log-weighted centroid -------------------------

ClusterHit HyCalCluster::reconstruct_pos(const ModuleCluster &cl) const
{
    const auto &center_mod = sys_.module(cl.center.index);

    // weights of the 3x3 neighbors relative to the total cluster energy
    LogCentroid acc(get_weight(cl.center.energy, cl.energy));
    int count = 0;
    for (auto &hit : cl.hits) {
        if (hit.index == cl.center.index) continue;
        if (count >= POS_RECON_HITS) break;

        double dx, dy;
        sys_.qdist(center_mod, sys_.module(hit.index), dx, dy);
        if (qdist_in_3x3(dx, dy)) {
            acc.add(dx, dy, get_weight(hit.energy, cl.energy));
            ++count;
        }
    }

    ClusterHit result;
    result.center_id = center_mod.id;
    result.energy    = cl.energy;
    result.time      = cl.center.time;
    result.nblocks   = static_cast<int>(cl.hits.size());
    result.flag      = cl.flag;
    result.linear_corr = 1.f;
    result.bias_corr = 1.f;
    result.leakage = cl.leakage;
    result.energy_square = cl.energy_square;

    acc.position(center_mod, result.x, result.y);
    result.npos = acc.n;

    // if available, update the weighted position with leakage correction
    if (cl.has_leakage_position) {
        result.x = cl.leakage_x;
        result.y = cl.leakage_y;
        result.npos = cl.leakage_npos;
    }

    if (config_.energy_bias_correction && config_.energy_bias) {
        result.bias_corr = config_.energy_bias->correction_factor(
            center_mod, result.x, result.y, cl.energy);
    }

    // non-linear correction as the last step
    if (config_.non_linear_corr) {
        // 1/linear_corr = E_rec/E_exp
        // = 1 + nl1*(E_rec-E_base)/1000 + nl2*((E_rec-E_base)/1000)^2
        const float nl1 = center_mod.cal_non_linear_1;
        const float nl2 = center_mod.cal_non_linear_2;
        const float base_energy = center_mod.cal_base_energy;
        const float delta_gev = (cl.energy - base_energy) / 1000.f;
        float non_linear_factor = 1.f / (1.f + nl1 * delta_gev
                                               + nl2 * delta_gev * delta_gev);
        if (cl.energy > 3800.f || non_linear_factor < 0.7f || non_linear_factor > 1.3f)
            non_linear_factor = 1.f;
        result.linear_corr *= non_linear_factor;
    }

    result.energy = cl.energy * result.bias_corr * result.linear_corr;
    result.energy_square *= result.bias_corr * result.linear_corr;

    return result;
}

float HyCalCluster::get_weight(float E, float E_total) const
{
    if (E_total <= 0.f) return 0.f;
    float w = config_.log_weight_thres + std::log(E / E_total);
    return (w > 0.f) ? w : 0.f;
}

void HyCalCluster::apply_leakage_correction(ModuleCluster &cl) const
{
    if (!config_.leakage_correction || test_bit(cl.flag, kLeakCorr)) return;
    if (cl.energy < config_.min_cluster_energy) return;
    if (static_cast<int>(cl.hits.size()) < config_.min_cluster_size) return;
    if (cl.hits.size() < 4) return;

    const auto &center_mod = sys_.module(cl.center.index);
    if (!test_bit(center_mod.flag, kLeakage) || center_mod.virtual_neighbors.empty())
        return;

    std::vector<LeakageHit> leaks;
    leaks.reserve(center_mod.virtual_neighbors.size());
    for (const auto &vn : center_mod.virtual_neighbors) {
        LeakageHit hit;
        hit.x = vn.x;
        hit.y = vn.y;
        hit.dx = vn.dx;
        hit.dy = vn.dy;
        hit.sector = center_mod.sector;
        hit.type = vn.type;
        if (vn.backing_module >= 0 && vn.backing_module < sys_.module_count())
            hit.sector = sys_.module(vn.backing_module).sector;
        leaks.push_back(hit);
    }
    if (leaks.empty()) return;

    LeakagePoint pos = reconstruct_leakage_position(cl, leaks, cl.energy);
    double est = eval_cluster_profile(pos, cl);
    if (!std::isfinite(est)) return;

    std::vector<float> previous(leaks.size(), 0.f);
    for (int iter = 0; iter < config_.leakage_iterations; ++iter) {
        for (size_t i = 0; i < leaks.size(); ++i)
            previous[i] = leaks[i].energy;

        float leakage = 0.f;
        for (auto &leak : leaks) {
            auto prof = get_pwo_profile_value_at(pos.x, pos.y, pos.energy,
                                                 leak.x, leak.y);
            leak.energy = 0.f;
            if (prof.frac > config_.least_leakage_fraction && prof.frac < 1.f)
                leak.energy = pos.energy * prof.frac;
            leakage += leak.energy;
        }

        if (!std::isfinite(leakage) || leakage <= 0.f) break;
        if (config_.max_leakage_fraction > 0.f &&
            leakage / cl.energy > config_.max_leakage_fraction) {
            // Reject the complete correction, including any earlier accepted
            // iteration, rather than leaving a partially divergent result.
            return;
        }

        LeakagePoint new_pos = reconstruct_leakage_position(cl, leaks,
                                                            cl.energy + leakage);
        double new_est = eval_cluster_profile(new_pos, cl);
        if (!std::isfinite(new_est) || new_est >= est) {
            for (size_t i = 0; i < leaks.size(); ++i)
                leaks[i].energy = previous[i];
            break;
        }

        const float relative_change = std::fabs(new_pos.energy - pos.energy) /
                                      std::max(pos.energy, 1.f);
        pos = new_pos;
        est = new_est;
        if (config_.leakage_convergence_rel > 0.f &&
            relative_change < config_.leakage_convergence_rel)
            break;
    }

    float leakage = 0.f;
    for (const auto &leak : leaks)
        leakage += leak.energy;
    if (leakage <= 0.f) return;

    cl.leakage = leakage;
    cl.energy += leakage;
    cl.energy_square += leakage;
    cl.has_leakage_position = true;
    cl.leakage_x = pos.x;
    cl.leakage_y = pos.y;
    cl.leakage_npos = pos.npos;
    set_bit(cl.flag, kLeakCorr);
}

HyCalCluster::LeakagePoint HyCalCluster::reconstruct_leakage_position(
    const ModuleCluster &cl,
    const std::vector<LeakageHit> &leaks,
    float total_energy) const
{
    const auto &center_mod = sys_.module(cl.center.index);
    LogCentroid acc(get_weight(cl.center.energy, total_energy));

    for (const auto &hit : cl.hits) {
        if (hit.index == cl.center.index) continue;

        double dx, dy;
        sys_.qdist(center_mod, sys_.module(hit.index), dx, dy);
        if (qdist_in_3x3(dx, dy))
            acc.add(dx, dy, get_weight(hit.energy, total_energy));
    }

    for (const auto &leak : leaks) {
        if (leak.energy > 0.f && qdist_in_3x3(leak.dx, leak.dy))
            acc.add(leak.dx, leak.dy, get_weight(leak.energy, total_energy));
    }

    LeakagePoint pos;
    acc.position(center_mod, pos.x, pos.y);
    pos.energy = total_energy;
    pos.npos = acc.n;
    return pos;
}

double HyCalCluster::eval_cluster_profile(const LeakagePoint &pos,
                                           const ModuleCluster &cl) const
{
    if (pos.energy <= 0.f) return std::numeric_limits<double>::infinity();

    const double sigma_E = sys_.EnergyResolution(pos.energy);
    double est = 0.;
    int count = 0;

    for (const auto &hit : cl.hits) {
        const auto &mod = sys_.module(hit.index);
        auto prof = get_pwo_profile_value_at(pos.x, pos.y, pos.energy,
                                             mod.x, mod.y);
        if (prof.frac < 0.01f) continue;

        const double diff = hit.energy - pos.energy * prof.frac;
        const double sigma2 = pos.energy * pos.energy * prof.err * prof.err +
                              sigma_E * sigma_E * prof.frac * prof.frac;
        if (sigma2 <= 0.) continue;

        est += std::abs(diff) / std::sqrt(sigma2);
        ++count;
    }

    return (count > 0) ? est / count : std::numeric_limits<double>::infinity();
}

// --- profile helpers --------------------------------------------------------

float HyCalCluster::get_profile_frac(const ModuleHit &center, const ModuleHit &hit) const
{
    const auto &m1 = sys_.module(center.index);
    const auto &m2 = sys_.module(hit.index);
    double dx, dy;
    sys_.qdist(m1, m2, dx, dy);
    float dist = std::sqrt(static_cast<float>(dx * dx + dy * dy));
    // center module contains ~78% of energy, scale up to estimate total
    return profile_->GetFraction(m1.type, dist, center.energy / 0.78f);
}

ProfileValue HyCalCluster::get_profile_value_at(float cx, float cy, float cE,
                                                 double mx, double my,
                                                 int msector,
                                                 ModuleType type) const
{
    int sid = sys_.get_sector_id(cx, cy);
    if (sid < 0 || sid >= static_cast<int>(Sector::Max))
        sid = msector;
    double dx, dy;
    sys_.qdist(cx, cy, sid, mx, my, msector, dx, dy);
    float dist = std::sqrt(static_cast<float>(dx * dx + dy * dy));
    return profile_->GetFractionValue(type, dist, cE);
}

ProfileValue HyCalCluster::get_pwo_profile_value_at(float cx, float cy, float cE,
                                                     double mx, double my) const
{
    // Leakage candidates are virtual continuations of the PbWO4 grid.  Use
    // the W cell pitch directly even when the reconstructed point moves past
    // the physical W boundary; sector-based qdist would otherwise introduce
    // PbGlass dimensions at the outer edge.
    const auto &w = sys_.sector_info(static_cast<int>(Sector::Center));
    if (w.msize_x <= 0. || w.msize_y <= 0.) return {};
    const double dx = (mx - cx) / w.msize_x;
    const double dy = (my - cy) / w.msize_y;
    const float dist = std::sqrt(static_cast<float>(dx * dx + dy * dy));
    return profile_->GetFractionValue(ModuleType::PbWO4, dist, cE);
}

float HyCalCluster::ProfileFractionAt(float cx, float cy, float cE,
                                      int module_index) const
{
    const auto &m = sys_.module(module_index);
    return get_profile_value_at(cx, cy, cE, m.x, m.y, m.sector, m.type).frac;
}

// Seed selection mirrors group_hits(); every other pulse within ±max_qdist on
// each axis is emitted, so a pulse may appear as a neighbour of more than one
// seed (the study is about the dt landscape, not cluster assignment).  Works
// on local scratch so the per-event state used by FormClusters() is untouched.
void HyCalCluster::CollectNeighborTiming(std::vector<SeedNeighborTiming> &out,
                                          double max_qdist) const
{
    out.clear();
    if (hits_.empty()) return;

    std::vector<std::vector<int>> mod_hits(sys_.module_count());
    for (int i = 0; i < static_cast<int>(hits_.size()); ++i)
        mod_hits[hits_[i].index].push_back(i);

    std::vector<int> order(hits_.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
              [&](int a, int b) { return hits_[a].energy > hits_[b].energy; });

    // A seed is "claimed" only against its own (module, time) so a single
    // physics shower with several pulses doesn't double-seed.  Use a small
    // claim window matching min_center_energy regions (5 ns is shorter than
    // any reasonable pulse separation) to keep the seed list compact.
    constexpr float SEED_CLAIM_WINDOW_NS = 5.f;

    std::vector<bool> seed_claimed(hits_.size(), false);

    for (int seed : order) {
        if (seed_claimed[seed]) continue;
        if (hits_[seed].energy < config_.min_center_energy) break;

        const float st = hits_[seed].time;
        const int   sm = hits_[seed].index;
        for (int k : mod_hits[sm])
            if (std::fabs(hits_[k].time - st) <= SEED_CLAIM_WINDOW_NS)
                seed_claimed[k] = true;

        const auto &smod = sys_.module(sm);
        for (int k = 0; k < static_cast<int>(hits_.size()); ++k) {
            if (k == seed) continue;
            const auto &nmod = sys_.module(hits_[k].index);
            double dx, dy;
            sys_.qdist(smod, nmod, dx, dy);
            if (std::fabs(dx) > max_qdist || std::fabs(dy) > max_qdist) continue;

            SeedNeighborTiming s;
            s.seed_module     = sm;
            s.neighbor_module = hits_[k].index;
            s.seed_time       = st;
            s.neighbor_time   = hits_[k].time;
            s.dt              = hits_[k].time - st;
            s.seed_energy     = hits_[seed].energy;
            s.neighbor_energy = hits_[k].energy;
            s.dx_q            = dx;
            s.dy_q            = dy;
            out.push_back(s);
        }
    }
}

} // namespace fdec
