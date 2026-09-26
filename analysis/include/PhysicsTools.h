#pragma once
//=============================================================================
// PhysicsTools.h — physics analysis tools for PRad2
//
// Provides kinematic calculations, energy loss corrections, and
// per-module energy histogram management with ROOT.
// Depends on prad2det (HyCalSystem) and ROOT (TH1F/TH2F).
//=============================================================================

#include "HyCalSystem.h"
#include <TH1F.h>
#include <TH2F.h>
#include <array>
#include <cmath>
#include <string>
#include <vector>
#include <memory>

class TH2Poly;
class TLorentzVector;

namespace analysis {

struct GEMHit {
    float x = 0.f;
    float y = 0.f;
    float z = 0.f;
    uint8_t det_id = 5; // 0-3 for GEM1-GEM4
};

struct HCHit {
    float x = 0.f;
    float y = 0.f;
    float z = 0.f;
    float energy = 0.f;
    uint16_t center_id = 0; // index of central block
    uint32_t flag = -1;
};

//data structure for storing reconstructed Moller events used for analysis
struct DataPoint
{
    float x;
    float y;
    float z;
    float E;

    DataPoint() {};
    DataPoint(float xi, float yi, float zi, float Ei) : x(xi), y(yi), z(zi), E(Ei) {};
};
typedef std::pair<DataPoint, DataPoint> MollerEvent;
typedef std::vector<MollerEvent> MollerData;

// Nominal PbWO4 module size (mm): size_y from hycal_map.json (size_x is 20.77).
// The fiducial ring bounds below are multiples of it in both x and y.
inline constexpr double kPbWO4Pitch = 20.75;

// Square-annulus fiducial cut: outside the inner square of half-width
// inner*pitch and inside the outer square of half-width outer*pitch.
// Frame-agnostic; the caller picks lab, HyCal or module-centre coordinates.
inline bool InHyCalRing(double x, double y, double inner, double outer,
                        double pitch = kPbWO4Pitch)
{
    return (std::fabs(x) > pitch * inner || std::fabs(y) > pitch * inner)
        && (std::fabs(x) < pitch * outer && std::fabs(y) < pitch * outer);
}

class PhysicsTools
{
public:
    explicit PhysicsTools(fdec::HyCalSystem &hycal);
    ~PhysicsTools();

    // --- per-module cluster energy histograms --------------------------------
    void FillModuleEnergy(int module_id, float energy);
    TH1F *GetModuleEnergyHist(int module_id) const;

    void FillEnergyVsModule(int module_id, float energy);
    TH2F *GetEnergyVsModuleHist() const { return h2_energy_module_.get(); }

    void FillEnergyVsTheta(float theta_deg, float energy);
    TH2F *GetEnergyVsThetaHist() const { return h2_energy_theta_.get(); }

    // must be called after all events are processed
    // and also need to call FillModuleEnergy for every event first
    void FillNeventsModuleMap() {
        for (int module_id = 1; module_id <= 1156; module_id++) {
            const auto *mod = hycal_.module_by_id(module_id + 1000);
            if (!mod || !mod->is_pwo4()) continue;
            TH1F *h = GetModuleEnergyHist(module_id + 1000);
            if (!h) continue;
            int count = h->GetEntries();
            h2_Nevents_moduleMap_->SetBinContent(mod->column + 1, 34 - mod->row, count);
        }
    }
    TH2F *GetNeventsModuleMapHist() const { return h2_Nevents_moduleMap_.get(); }

    // TH2Poly with one rectangular bin per PbWO4 module, axes ±half_range (mm),
    // created with `new` in the current directory (caller owns it).
    // bin_by_index is resized to module_count(): the TH2Poly bin of each
    // PbWO4 module by module index, -1 for every other module.
    static TH2Poly *MakeModuleMap(const fdec::HyCalSystem &hycal, const char *name,
                                  const char *title, double half_range,
                                  std::vector<int> &bin_by_index);

    // --- peak / resolution analysis ------------------------------------------
    // Returns {peak, sigma, chi2} from Gaussian fit.
    std::array<float, 3> FitPeakResolution(int module_id) const;
    static std::array<double, 5> fitGaus(TH1F *h, float expectPeak = 0.f,
                                        bool withError = false);
    static std::array<double, 5> fitCrystalBall(TH1F *h, float expectPeak = 0.f,
                                               float alpha = 1.5f, float n = 5.0f,
                                               bool withError = false);
    // Returns {mean, sigma, chi2/ndf, mean_error, sigma_error}.
    // The last two values are zero unless withError is true.
    static std::array<double, 5> fitPeak(TH1F *h, float expectPeak = 0.f, bool withError = false,
                                        bool useCrystalBall = false,
                                        float alpha = 0.5f, float n = 5.0f);

    // --- kinematics ----------------------------------------------------------
    static constexpr float kProtonMass   = 938.272f;     // MeV
    static constexpr float kElectronMass = 0.51099895f;  // MeV

    // Expected energy for elastic e-p or e-e scattering.
    //   theta: scattering angle in degrees
    //   Ebeam: beam energy in MeV
    //   type:  "ep" or "ee"
    static float ExpectedEnergy(float theta_deg, float Ebeam, const std::string &type);

    // Four-momentum of a particle of mass m (MeV) and energy E (MeV) emitted
    // from the target (origin) towards the hit at (x, y, z).  Returns false and
    // zeroes p4 when E < m or the hit is at the origin.
    static bool HitP4(float x, float y, float z, float E, float m, TLorentzVector &p4);

    // Energy loss correction for electron passing through target + windows.
    //   theta: scattering angle in degrees
    //   E:     measured energy in MeV
    static float EnergyLoss(float theta_deg, float E);

    // elastic e-e kinematic check for Moller event selection: energy sum
    // within 5 sigma of EBeam and each energy within 3.5 sigma of its
    // expected value, sigma = resolution * E / sqrt(E in GeV)
    static bool isMoller_kinematic(float theta_deg1, float energy1, float theta_deg2, float energy2, float EBeam, float resolution);

    //physics analysis helpers

    // Get the center of the Moller distribution in x-y space
    // enter two moller events, find the intersection of 2 lines
    // output the x-y coordinates of the center for every 2 moller events
    static std::array<float, 2> GetMollerCenter(const MollerEvent &event1,
                                                const MollerEvent &event2);

    static float GetMollerZdistance(const MollerEvent &event, float Ebeam);

    //Get azimuthal angle difference(should be around 180 degrees) of the Moller event
    static float GetMollerPhiDiff(const MollerEvent &event1);

    // The two hits are back to back in phi within max_dev_deg.
    static bool isBackToBack(const MollerEvent &event, float max_dev_deg)
    {
        return std::fabs(GetMollerPhiDiff(event)) < max_dev_deg;
    }

    static float GetPhiAngle(float x, float y);

    // Polar angle (degrees) of (x, y, z) seen from the target (origin).
    static float GetThetaAngle(float x, float y, float z);

    void FillMollerPhiDiff(float phi_diff) { if (moller_phi_diff_) moller_phi_diff_->Fill(phi_diff); }
    void FillMollerXY(float x, float y) { if (moller_x_) moller_x_->Fill(x); if (moller_y_) moller_y_->Fill(y); }
    void FillMollerZ(float z) { if (moller_z_) moller_z_->Fill(z); }

    TH1F *GetMollerPhiDiffHist() const { return moller_phi_diff_.get(); };
    TH1F *GetMollerXHist() const { return moller_x_.get(); };
    TH1F *GetMollerYHist() const { return moller_y_.get(); };
    TH1F *GetMollerZHist() const { return moller_z_.get(); };

private:
    fdec::HyCalSystem &hycal_;
    std::vector<std::unique_ptr<TH1F>> module_hists_;  // one per module
    std::unique_ptr<TH2F> h2_energy_module_;
    std::unique_ptr<TH2F> h2_energy_theta_;
    std::unique_ptr<TH2F> h2_Nevents_moduleMap_;
    std::unique_ptr<TH1F> moller_phi_diff_;
    std::unique_ptr<TH1F> moller_x_;
    std::unique_ptr<TH1F> moller_y_;
    std::unique_ptr<TH1F> moller_z_;
};

} // namespace analysis
