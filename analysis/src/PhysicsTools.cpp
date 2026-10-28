//=============================================================================
// PhysicsTools.cpp — physics analysis tools
//=============================================================================

#include "PhysicsTools.h"
#include <TF1.h>
#include <TH2Poly.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <algorithm>
#include <cmath>

namespace analysis {

static constexpr float DEG2RAD = 3.14159265f / 180.f;

PhysicsTools::PhysicsTools(fdec::HyCalSystem &hycal)
    : hycal_(hycal)
{
    int nmod = hycal_.module_count();
    module_hists_.resize(nmod);
    for (int i = 0; i < nmod; ++i) {
        auto &mod = hycal_.module(i);
        std::string name = "h_" + mod.name;
        std::string title = mod.name + " cluster energy;Energy (MeV);Counts";
        module_hists_[i] = std::make_unique<TH1F>(name.c_str(), title.c_str(), 500, 0, 5000);
    }
    h2_energy_module_ = std::make_unique<TH2F>(
        "h2_energy_module", "Energy vs Module;Module Index;Energy (MeV)",
        nmod, 0, nmod, 2000, 0, 4000);
    h2_energy_theta_ = std::make_unique<TH2F>(
        "h2_energy_theta", "Energy vs Theta;Theta (deg);Energy (MeV)",
        160, 0, 8, 4000, 0, 4000);
    h2_Nevents_moduleMap_ = std::make_unique<TH2F>(
        "h2_Nevents_moduleMap", "Number of Events per Module;Column;Row",
        34, 0.5, 34.5, 34, -34.5, -0.5);

    moller_phi_diff_ = std::make_unique<TH1F>(
        "h_moller_phi_diff", "Moller Phi Difference;Phi_{e1} - Phi_{e2} (deg);Counts",
        40, -20, 20);
    moller_x_ = std::make_unique<TH1F>(
        "h_moller_x", "Moller Center X Position (HyCal);X (mm);Counts",
        100, -10, 10);
    moller_y_ = std::make_unique<TH1F>(
        "h_moller_y", "Moller Center Y Position (HyCal);Y (mm);Counts",
        100, -10, 10);
    moller_z_ = std::make_unique<TH1F>(
        "h_moller_z", "Moller Z Position (HyCal);Z (mm);Counts",
        1000, 5000, 8000);
}

PhysicsTools::~PhysicsTools() = default;

namespace {

// Crystal Ball: p[0]=amp, p[1]=mean, p[2]=sigma, p[3]=alpha, p[4]=n
double crystalBallFunc(double *x, double *p)
{
    double amp   = p[0];
    double mu    = p[1];
    double sigma = p[2];
    double alpha = p[3];
    double n     = p[4];
    double t = (x[0] - mu) / sigma;
    if (t > -std::abs(alpha)) {
        return amp * std::exp(-0.5 * t * t);
    } else {
        double a = std::pow(n / std::abs(alpha), n) * std::exp(-0.5 * alpha * alpha);
        double b = n / std::abs(alpha) - std::abs(alpha);
        return amp * a * std::pow(b - t, -n);
    }
}

// Fit window of fitGaus / fitCrystalBall.
struct PeakWindow {
    double height = 0.;   // content of the peak bin
    double peak0  = 0.;   // centre of the peak bin
    double sigma0 = 0.;   // expected width, 3% / sqrt(E in GeV)
    double lo = 0., hi = 0.;
};

// Peak: the largest local maximum with its bin in [loMul, hiMul] * expectPeak,
// else the global maximum.  Window: grown from the peak while bins exceed
// threshFrac of its height, widened to at least 2.5 sigma0 and shrunk to at
// most maxSpan sigma0, from bin centre to bin centre.  Returns false when
// the histogram is not fittable (< 100 entries, < 4 bins, empty peak or
// window).  Sets Poisson bin errors (Sumw2) for the chi-square fit.
bool findPeakWindow(TH1F *h, float expectPeak, double loMul, double hiMul,
                    double threshFrac, double maxSpan, PeakWindow &w)
{
    if (!h || h->GetEntries() < 100) return false;

    const int nBins = h->GetNbinsX();
    if (nBins < 4) return false;

    int peakBin = -1;
    double peakHeight = 0.;

    if (std::isfinite(expectPeak) && expectPeak > 0.) {
        int firstBin = h->GetXaxis()->FindFixBin(loMul * expectPeak);
        int lastBin  = h->GetXaxis()->FindFixBin(hiMul * expectPeak);
        firstBin = std::max(1, firstBin);
        lastBin  = std::min(nBins, lastBin);

        for (int bin = firstBin; bin <= lastBin; ++bin) {
            const double content = h->GetBinContent(bin);
            const double left = (bin > 1) ? h->GetBinContent(bin - 1) : content;
            const double right = (bin < nBins) ? h->GetBinContent(bin + 1) : content;
            if (content > 0. && content >= left && content >= right && content > peakHeight) {
                peakBin = bin;
                peakHeight = content;
            }
        }
    }

    if (peakBin < 0) {
        peakBin = h->GetMaximumBin();
        peakHeight = h->GetBinContent(peakBin);
    }
    if (peakHeight <= 0.) return false;

    const double threshold = threshFrac * peakHeight;
    const double peak0 = h->GetBinCenter(peakBin);
    const double sigma0 = 0.03 * sqrt(peak0 * 1000.);
    int leftBin = peakBin;
    int rightBin = peakBin;
    while (leftBin > 1 && h->GetBinContent(leftBin) > threshold) --leftBin;
    while (rightBin < nBins && h->GetBinContent(rightBin) > threshold) ++rightBin;
    while ((h->GetBinCenter(rightBin) - h->GetBinCenter(leftBin)) < 2.5 * sigma0
           && (leftBin > 1 || rightBin < nBins)) {
        if (leftBin > 1) --leftBin;
        if (rightBin < nBins) ++rightBin;
    }
    while ((h->GetBinCenter(rightBin) - h->GetBinCenter(leftBin)) > maxSpan * sigma0
           && (leftBin > 1 || rightBin < nBins)) {
        ++leftBin;
        --rightBin;
    }

    const double lo = h->GetBinCenter(leftBin);
    const double hi = h->GetBinCenter(rightBin);
    if (!(hi > lo) || !std::isfinite(sigma0) || sigma0 <= 0.) return false;

    // ROOT's chi-square fit uses the histogram bin errors. Sumw2 initializes
    // Poisson statistical errors for an unweighted histogram and preserves them
    // correctly if the histogram is filled again later.
    if (h->GetSumw2N() == 0) h->Sumw2();

    w = {peakHeight, peak0, sigma0, lo, hi};
    return true;
}

// Fit f over its range and return {mean, sigma, chi2/ndf, mean_error,
// sigma_error} from its parameters 1 and 2, all zero when the fit failed.
// The fitted TF1 stays attached to the histogram, so it is drawn with it and
// persisted when the histogram is written to a ROOT file.
std::array<double, 5> fitPeakModel(TH1F *h, TF1 &f, bool withError)
{
    const int fitStatus = h->Fit(&f, "RQ");
    if (fitStatus != 0) return {0., 0., 0., 0., 0.};

    const double mean = f.GetParameter(1);
    const double sigma = std::abs(f.GetParameter(2));
    if (!std::isfinite(mean) || !std::isfinite(sigma) || sigma <= 0.)
        return {0., 0., 0., 0., 0.};

    const double chi2 = (f.GetNDF() > 0) ? f.GetChisquare() / f.GetNDF() : 0.;
    return {mean, sigma, chi2,
            withError ? f.GetParError(1) : 0.,
            withError ? f.GetParError(2) : 0.};
}

} // anonymous namespace

void PhysicsTools::FillModuleEnergy(int module_id, float energy)
{   
    if (module_id >= 0){
        int module_index = hycal_.id_to_index(module_id);
        if (module_index >= 0 && module_index < (int)module_hists_.size())
            module_hists_[module_index]->Fill(energy);
    }
}

TH1F *PhysicsTools::GetModuleEnergyHist(int module_id) const
{
    int module_index = hycal_.id_to_index(module_id);
    if (module_index >= 0 && module_index < (int)module_hists_.size())
        return module_hists_[module_index].get();
    return nullptr;
}

void PhysicsTools::FillEnergyVsModule(int module_id, float energy)
{
    int module_index = hycal_.id_to_index(module_id);
    if (module_index >= 0 && module_index < (int)module_hists_.size())
        h2_energy_module_->Fill(module_index, energy);
}

void PhysicsTools::FillEnergyVsTheta(float theta_deg, float energy)
{
    if (h2_energy_theta_)
        h2_energy_theta_->Fill(theta_deg, energy);
}

TH2Poly *PhysicsTools::MakeModuleMap(const fdec::HyCalSystem &hycal, const char *name,
                                     const char *title, double half_range,
                                     std::vector<int> &bin_by_index)
{
    auto *poly = new TH2Poly(name, title, -half_range, half_range, -half_range, half_range);
    bin_by_index.assign(hycal.module_count(), -1);
    for (int m = 0; m < hycal.module_count(); ++m) {
        const auto &mod = hycal.module(m);
        if (!mod.is_pwo4()) continue;
        bin_by_index[m] = poly->AddBin(mod.x - 0.5 * mod.size_x, mod.y - 0.5 * mod.size_y,
                                       mod.x + 0.5 * mod.size_x, mod.y + 0.5 * mod.size_y);
    }
    return poly;
}

std::array<float, 3> PhysicsTools::FitPeakResolution(int module_id) const
{
    int module_index = hycal_.id_to_index(module_id);
    if (module_index < 0 || module_index >= (int)module_hists_.size())
        return {0.f, 0.f, 0.f};

    TH1F *h = module_hists_[module_index].get();
    if (!h || h->GetEntries() < 1) return {0.f, 0.f, 100.f};

    const double resolution = 0.035; // 3.5% / sqrt(E/1000) energy resolution

    auto estimateSigma = [&](double E) -> double {
        return (E > 0.) ? E * resolution / std::sqrt(E / 1000.) : 1.;
    };

    // Step 1: find histogram maximum as initial center
    double center = h->GetBinCenter(h->GetMaximumBin());

    // Step 2: weighted mean within [center ± 3σ]
    double sigma = estimateSigma(center);
    double sumW = 0., sumWx = 0.;
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
        double x = h->GetBinCenter(b);
        if (x < center - 3.*sigma || x > center + 3.*sigma) continue;
        double w = h->GetBinContent(b);
        sumW  += w;
        sumWx += w * x;
    }
    double mean = (sumW > 0.) ? sumWx / sumW : center;

    // Step 3: first Gaussian fit within [mean ± 2σ], σ re-estimated from new mean
    sigma = estimateSigma(mean);
    {
        TF1 g1("_fpk_g1_", "gaus", mean - 2.*sigma, mean + 2.*sigma);
        g1.SetParameters(h->GetMaximum(), mean, sigma);
        h->Fit(&g1, "RQ0");
        mean  = g1.GetParameter(1);
        // keep sigma from resolution for next iteration boundary
    }

    // Step 4: final Gaussian fit within [mean ± 1σ], σ re-estimated from new mean
    sigma = estimateSigma(mean);
    TF1 g2("_fpk_g2_", "gaus", mean - sigma, mean + sigma);
    g2.SetParameters(h->GetMaximum(), mean, sigma);
    h->Fit(&g2, "RQ0");
    mean  = g2.GetParameter(1);
    sigma = std::abs(g2.GetParameter(2));
    double chi2 = (g2.GetNDF() > 0) ? g2.GetChisquare() / g2.GetNDF() : 0.;

    return {static_cast<float>(mean), static_cast<float>(sigma), static_cast<float>(chi2)};
}

float PhysicsTools::ExpectedEnergy(float theta_deg, float Ebeam, const std::string &type)
{
    float theta = theta_deg * DEG2RAD;
    float cos_t = std::cos(theta);
    float sin_t = std::sin(theta);

    if (type == "ep") {
        // elastic e-p: E' = E * M / (M + E*(1 - cos_t))
        // where M = proton mass
        float expectE = Ebeam * kProtonMass / (kProtonMass + Ebeam * (1.f - cos_t));
        float eloss = EnergyLoss(theta_deg, expectE);
        return expectE - eloss;
    }
    if (type == "ee") {
        // Moller scattering: exact lab-frame formula from 4-momentum conservation
        // E' = m * [(gamma+1) + (gamma-1)*cos^2(theta)] / [(gamma+1) - (gamma-1)*cos^2(theta)]
        float gamma = Ebeam / kElectronMass;
        float num = (gamma + 1.f) + (gamma - 1.f) * cos_t * cos_t;
        float den = (gamma + 1.f) - (gamma - 1.f) * cos_t * cos_t;
        if (den <= 0) return 0.f;
        float expectE = kElectronMass * num / den;
        float eloss = EnergyLoss(theta_deg, expectE);
        return expectE - eloss;
    }
    return 0.f;
}

float PhysicsTools::EnergyLoss(float theta_deg, float E)
{
    // simplified energy loss through target materials
    // path lengths scale as 1/cos(theta) for small angles
    float theta = theta_deg * DEG2RAD;
    float sec = (std::cos(theta) > 0.01f) ? (1.f / std::cos(theta)) : 100.f;

    // material thicknesses (mm) and dE/dx (MeV/mm) — approximate values
    // aluminum window: 0.5 mm, dE/dx ~ 1.6 MeV/mm, vacuum window
    // GEM foils: ~0.05 mm effective per GEM, dE/dx ~ 2.0 MeV/mm
    // GEM win Al foils: ~0.06 mm effective per GEM, dE/dx ~ 1.6 MeV/mm
    // kapton window: ~0.24 mm per GEM, dE/dx ~ 1.8 MeV/mm,
    float eloss = 0.f;
    eloss += 0.500f * 1.6f * sec;  // Al window
    eloss += 0.120f * 1.6f * sec;  // GEM win Al foils (2 GEMs)
    eloss += 0.100f * 2.0f * sec;  // GEM foils (2 GEMs)
    eloss += 0.480f * 1.8f * sec;  // kapton cover

    return eloss;  // total energy loss in MeV
}

bool PhysicsTools::HitP4(float x, float y, float z, float E, float m, TLorentzVector &p4)
{
    const float norm = std::sqrt(x * x + y * y + z * z);
    if (E < m || norm <= 0.f) {
        p4.SetXYZT(0., 0., 0., 0.);
        return false;
    }
    const float p = std::sqrt(E * E - m * m);
    p4.SetXYZT(p * (x / norm), p * (y / norm), p * (z / norm), E);
    return true;
}

bool PhysicsTools::isMoller_kinematic(float theta_deg1, float energy1, float theta_deg2, float energy2, float EBeam, float resolution)
{
    float expectE1 = ExpectedEnergy(theta_deg1, EBeam, "ee");
    float expectE2 = ExpectedEnergy(theta_deg2, EBeam, "ee");

    return fabs(energy1 + energy2 - EBeam) < 5.f * resolution * EBeam / sqrt(EBeam/1000.f)
        && fabs(energy1 - expectE1) < 3.5f * expectE1 * resolution / sqrt(expectE1/1000.f)
        && fabs(energy2 - expectE2) < 3.5f * expectE2 * resolution / sqrt(expectE2/1000.f);
}

std::array<float, 2> PhysicsTools::GetMollerCenter(const MollerEvent &event1,
                                                   const MollerEvent &event2)
{
    float x1[2], y1[2];
    float x2[2], y2[2];

    x1[0] = event1.first.x; y1[0] = event1.first.y;
    x1[1] = event1.second.x; y1[1] = event1.second.y;
    x2[0] = event2.first.x; y2[0] = event2.first.y;
    x2[1] = event2.second.x; y2[1] = event2.second.y;

    //two lines: y = ax + b, y = cx + d
    float dx1 = x1[0] - x1[1];
    float dx2 = x2[0] - x2[1];
    if (std::abs(dx1) < 1e-6f || std::abs(dx2) < 1e-6f)
        return {0.f, 0.f};  // vertical line — degenerate

    float a = (y1[0] - y1[1]) / dx1;
    float b = y1[0] - a * x1[0];
    float c = (y2[0] - y2[1]) / dx2;
    float d = y2[0] - c * x2[0];

    if (std::abs(a - c) < 1e-6f)
        return {0.f, 0.f};  // parallel lines — no intersection

    float x_cross = (d - b) / (a - c);
    float y_cross = a * x_cross + b;

    return {x_cross, y_cross};

}

float PhysicsTools::GetMollerZdistance(const MollerEvent &event, float Ebeam)
{
    float R1 = sqrt(event.first.x*event.first.x + event.first.y*event.first.y);
    float R2 = sqrt(event.second.x*event.second.x + event.second.y*event.second.y);
    float z = sqrt( (Ebeam + kElectronMass) * R1 * R2 / (2.*kElectronMass) );
    return z;
}

float PhysicsTools::GetMollerPhiDiff(const MollerEvent &event1)
{
    float x1 = event1.first.x, y1 = event1.first.y;
    float x2 = event1.second.x, y2 = event1.second.y;
    float phi1 = GetPhiAngle(x1, y1);
    float phi2 = GetPhiAngle(x2, y2);
    float phi_diff = fabs(phi1 - phi2) - 180.f; // Expecting back-to-back, so difference should be around 180 degrees
    return phi_diff;
}

float PhysicsTools::GetPhiAngle(float x, float y)
{
    // atan2 handles all quadrants and x==0 correctly
    float phi = std::atan2(y, x) * 180.f / static_cast<float>(TMath::Pi());
    if (phi < 0) phi += 360.f;
    return phi;
}

float PhysicsTools::GetThetaAngle(float x, float y, float z)
{
    return std::atan2(std::sqrt(x * x + y * y), z) * 180.0 / M_PI;
}

std::array<double, 5> PhysicsTools::fitPeak(TH1F *h, float expectPeak, bool withError,
                                            bool useCrystalBall,
                                            float alpha, float n)
{
    if (useCrystalBall) {
        return fitCrystalBall(h, expectPeak, alpha, n, withError);
    }
    return fitGaus(h, expectPeak, withError);
}

// Fit a peak near expectPeak with a Gaussian and return
// {mean, sigma, chi2/ndf, mean_error, sigma_error}.
// Returns all zeroes if the histogram or fit is invalid.
std::array<double, 5> PhysicsTools::fitGaus(TH1F *h, float expectPeak, bool withError)
{
    PeakWindow w;
    if (!findPeakWindow(h, expectPeak, 0.8, 1.2, 0.4, 4.5, w)) return {0., 0., 0., 0., 0.};

    TF1 gaus("_fg_", "gaus", w.lo, w.hi);
    gaus.SetParameters(w.height, w.peak0, w.sigma0);
    return fitPeakModel(h, gaus, withError);
}

// Fit a peak near expectPeak with a Crystal Ball and return
// {mean, sigma, chi2/ndf, mean_error, sigma_error}.
// Returns all zeroes if the histogram or fit is invalid.
std::array<double, 5> PhysicsTools::fitCrystalBall(TH1F *h, float expectPeak,
                                                  float alpha, float n, bool withError)
{
    PeakWindow w;
    if (!findPeakWindow(h, expectPeak, 0.7, 1.3, 0.05, 6.0, w)) return {0., 0., 0., 0., 0.};

    TF1 cb("_fcb_", crystalBallFunc, w.lo, w.hi, 5);
    cb.SetParName(0, "amp");
    cb.SetParName(1, "mean");
    cb.SetParName(2, "sigma");
    cb.SetParName(3, "alpha");
    cb.SetParName(4, "n");
    cb.SetParameters(w.height, w.peak0, w.sigma0, alpha, n);
    cb.SetParLimits(0, 0.0, 5.0 * w.height);
    cb.SetParLimits(1, w.lo, w.hi);
    cb.SetParLimits(2, 1e-6, std::max(w.hi - w.lo, 1e-3));
    cb.SetParLimits(3, 0.5, 4.0);
    cb.SetParLimits(4, 1.01, 20.0);
    return fitPeakModel(h, cb, withError);
}

}
