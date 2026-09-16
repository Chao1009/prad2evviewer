// Crystal Ball: p[0]=amp, p[1]=mean, p[2]=sigma, p[3]=alpha, p[4]=n
static double crystalBallFunc(double *x, double *p)
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
// Gaussian: p[0]=amp, p[1]=mean, p[2]=sigma
static double gaussianFunc(double *x, double *p)
{
    double amp   = p[0];
    double mu    = p[1];
    double sigma = p[2];
    double t = (x[0] - mu) / sigma;
    return amp * std::exp(-0.5 * t * t);
}

// Fit a peak near expectPeak with a Gaussian and return {mean, sigma, chi2/ndf}.
// Returns {0,0,0} if the histogram or fit is invalid.
std::array<double, 3> PhysicsTools::fitGaus(TH1F *h, float expectPeak)
{
    if (!h || h->GetEntries() < 100) return {0., 0., 0.};

    const int nBins = h->GetNbinsX();
    if (nBins < 4) return {0., 0., 0.};

    int peakBin = -1;
    double peakHeight = 0.;

    // Prefer the largest local maximum within +/-20% of the expected peak.
    if (std::isfinite(expectPeak) && expectPeak > 0.) {
        int firstBin = h->GetXaxis()->FindFixBin(0.8 * expectPeak);
        int lastBin  = h->GetXaxis()->FindFixBin(1.2 * expectPeak);
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

    // Fall back to the global maximum when no peak is found near expectPeak.
    if (peakBin < 0) {
        peakBin = h->GetMaximumBin();
        peakHeight = h->GetBinContent(peakBin);
    }
    if (peakHeight <= 0.) return {0., 0., 0.};

    const double threshold = 0.4 * peakHeight;
    int leftBin = peakBin;
    int rightBin = peakBin;
    while (leftBin > 1 && h->GetBinContent(leftBin) > threshold) --leftBin;
    while (rightBin < nBins && h->GetBinContent(rightBin) > threshold) ++rightBin;
    // Widen a too-narrow threshold-crossing range to a minimum 50 MeV span
    // so low-statistics spectra still get a usable fit window.
    while ((h->GetBinCenter(rightBin) - h->GetBinCenter(leftBin)) < 50.
           && (leftBin > 1 || rightBin < nBins)) {
        if (leftBin > 1) --leftBin;
        if (rightBin < nBins) ++rightBin;
    }
    if (rightBin - leftBin + 1 < 4) return {0., 0., 0.};

    const double lo = h->GetBinCenter(leftBin);
    const double hi = h->GetBinCenter(rightBin);
    const double peak0 = h->GetBinCenter(peakBin);
    const double sigma0 = (hi - lo) / (2. * std::sqrt(-2. * std::log(0.4)));
    if (!(hi > lo) || !std::isfinite(sigma0) || sigma0 <= 0.) return {0., 0., 0.};

    // ROOT's chi-square fit uses the histogram bin errors. Sumw2 initializes
    // Poisson statistical errors for an unweighted histogram and preserves them
    // correctly if the histogram is filled again later.
    if (h->GetSumw2N() == 0) h->Sumw2();

    TF1 gaus("_fg_", "gaus", lo, hi);
    gaus.SetParameters(peakHeight, peak0, sigma0);
    const int fitStatus = h->Fit(&gaus, "RQN");
    if (fitStatus != 0) return {0., 0., 0.};

    const double mean = gaus.GetParameter(1);
    const double sigma = std::abs(gaus.GetParameter(2));
    if (!std::isfinite(mean) || !std::isfinite(sigma) || sigma <= 0.)
        return {0., 0., 0.};

    double chi2 = (gaus.GetNDF() > 0) ? gaus.GetChisquare() / gaus.GetNDF() : 0.;
    return {mean, sigma, chi2};
}

// Fit a peak near expectPeak with a Crystal Ball and return {mean, sigma, chi2/ndf}.
// Returns {0,0,0} if the histogram or fit is invalid.
std::array<double, 3> PhysicsTools::fitCrystalBall(TH1F *h, float expectPeak,
                                                  float alpha, float n)
{
    if (!h || h->GetEntries() < 100) return {0., 0., 0.};

    const int nBins = h->GetNbinsX();
    if (nBins < 4) return {0., 0., 0.};

    int peakBin = -1;
    double peakHeight = 0.;

    if (std::isfinite(expectPeak) && expectPeak > 0.) {
        int firstBin = h->GetXaxis()->FindFixBin(0.7 * expectPeak);
        int lastBin  = h->GetXaxis()->FindFixBin(1.3 * expectPeak);
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
    if (peakHeight <= 0.) return {0., 0., 0.};

    const double threshold = 0.1 * peakHeight;
    int leftBin = peakBin;
    int rightBin = peakBin;
    while (leftBin > 1 && h->GetBinContent(leftBin) > threshold) --leftBin;
    while (rightBin < nBins && h->GetBinContent(rightBin) > threshold) ++rightBin;
    while ((h->GetBinCenter(rightBin) - h->GetBinCenter(leftBin)) < 100.
           && (leftBin > 1 || rightBin < nBins)) {
        if (leftBin > 1) --leftBin;
        if (rightBin < nBins) ++rightBin;
    }
    if (rightBin - leftBin + 1 < 4) return {0., 0., 0.};

    const double lo = h->GetBinCenter(leftBin);
    const double hi = h->GetBinCenter(rightBin);
    const double peak0 = h->GetBinCenter(peakBin);
    const double sigma0 = (hi - lo) / (2. * std::sqrt(-2. * std::log(0.4)));
    if (!(hi > lo) || !std::isfinite(sigma0) || sigma0 <= 0.) return {0., 0., 0.};

    if (h->GetSumw2N() == 0) h->Sumw2();

    TF1 cb("_fcb_", crystalBallFunc, lo, hi, 5);
    cb.SetParName(0, "amp");
    cb.SetParName(1, "mean");
    cb.SetParName(2, "sigma");
    cb.SetParName(3, "alpha");
    cb.SetParName(4, "n");
    cb.SetParameters(peakHeight, peak0, sigma0, alpha, n);
    cb.SetParLimits(0, 0.0, std::numeric_limits<double>::max());
    cb.SetParLimits(1, lo, hi);
    cb.SetParLimits(2, 1e-6, std::max(hi - lo, 1e-3));
    cb.SetParLimits(3, 0.5, 30.0);
    cb.SetParLimits(4, 1.1, 100.0);

    const int fitStatus = h->Fit(&cb, "RQN");
    if (fitStatus != 0) return {0., 0., 0.};

    const double mean = cb.GetParameter(1);
    const double sigma = std::abs(cb.GetParameter(2));
    if (!std::isfinite(mean) || !std::isfinite(sigma) || sigma <= 0.)
        return {0., 0., 0.};

    const double chi2 = (cb.GetNDF() > 0) ? cb.GetChisquare() / cb.GetNDF() : 0.;
    return {mean, sigma, chi2};
}

void peak_fitting_compare() {
    // Implement the peak fitting comparison logic here.

    // a plot to draw the comparison between Gaussian and Crystal Ball curves
    // a gaus curve and several crystal ball curves with same amplitude and sigma 
    // but with different parameters of alpha and n
    TCanvas *curves = new TCanvas("curves", "Peak Fitting Comparison", 800, 600);

    const double amp = 1000.;
    const double mu = 2239.;
    const double sigma = 0.03 * sqrt(mu * 1000.);

    TF1 *gaus = new TF1("gaus_curve", gaussianFunc, mu - 6 * sigma, mu + 6 * sigma, 3);
    gaus->SetParameters(amp, mu, sigma);
    gaus->SetLineColor(kBlack);
    gaus->SetLineWidth(2);
    gaus->SetTitle("Gaussian vs Crystal Ball");
    gaus->GetXaxis()->SetTitle("Energy");
    gaus->GetYaxis()->SetTitle("Yield");

    std::vector<std::pair<std::string, std::pair<double, double>>> cb_configs = {
        {"CB #alpha=0.5, n=5", {0.5, 5.0}},
        {"CB #alpha=1.0, n=5", {1.0, 5.0}},
        {"CB #alpha=1.5, n=5", {1.5, 5.0}},
        {"CB #alpha=3.0, n=5", {3.0, 5.0}}
    };

    std::vector<TF1*> cbs;
    for (size_t i = 0; i < cb_configs.size(); ++i) {
        auto alpha = cb_configs[i].second.first;
        auto n = cb_configs[i].second.second;
        auto *cb = new TF1(Form("cb_%zu", i), crystalBallFunc, mu - 6 * sigma, mu + 6 * sigma, 5);
        cb->SetParameters(amp, mu, sigma, alpha, n);
        cb->SetLineColor(i+1);
        cb->SetLineWidth(2);
        cb->SetLineStyle(2);
        cbs.push_back(cb);
    }

    gaus->Draw();
    for (auto *cb : cbs) {
        cb->Draw("same");
    }

    gaus->GetXaxis()->SetRangeUser(mu - 6. * sigma, mu + 6. * sigma);
    gaus->GetYaxis()->SetRangeUser(0., 1.1 * amp);

    auto *leg = new TLegend(0.6, 0.55, 0.9, 0.9);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->AddEntry(gaus, "Gaussian", "l");
    for (size_t i = 0; i < cbs.size(); ++i) {
        leg->AddEntry(cbs[i], cb_configs[i].first.c_str(), "l");
    }
    leg->Draw();

    curves->Update();

    // Second canvas: compare Gaussian fit and Crystal Ball fit on the same histogram
    TCanvas *fit_compare = new TCanvas("fit_compare", "Gaussian vs Crystal Ball Fit", 900, 600);
    fit_compare->cd();

    TFile *file = TFile::Open("energy_spectra.root");
    TH1F *h1_Espec = nullptr;
    if (file && !file->IsZombie()) {
        h1_Espec = dynamic_cast<TH1F*>(file->Get("h1_Espec"));
    }
    if (!h1_Espec) {
        h1_Espec = new TH1F("h1_Espec_demo", "Demo energy spectrum;Energy;Counts",
                            120, 0.0, 5000.0);
        for (int i = 0; i < 200000; ++i) {
            double x = gRandom->Gaus(2239.0, 80.0);
            if (x > 0.0) h1_Espec->Fill(x);
        }
    }

    h1_Espec->SetLineColor(kGray + 2);
    h1_Espec->Draw("hist");

    double fit_min = h1_Espec->GetXaxis()->GetXmin();
    double fit_max = h1_Espec->GetXaxis()->GetXmax();

    TF1 *gaus_fit = new TF1("gaus_fit", "gaus", fit_min, fit_max);
    gaus_fit->SetLineColor(kBlue);
    gaus_fit->SetLineWidth(2);
    h1_Espec->Fit(gaus_fit, "RQN");

    TF1 *cb_fit = new TF1("cb_fit", crystalBallFunc, fit_min, fit_max, 5);
    cb_fit->SetParameters(gaus_fit->GetParameter(0),
                          gaus_fit->GetParameter(1),
                          gaus_fit->GetParameter(2),
                          1.5, 3.0);
    cb_fit->SetParLimits(0, 0.0, 1.0e30);
    cb_fit->SetParLimits(1, h1_Espec->GetMean() - 3.0 * h1_Espec->GetRMS(),
                         h1_Espec->GetMean() + 3.0 * h1_Espec->GetRMS());
    cb_fit->SetParLimits(2, 1.0e-6, 1.0e4);
    cb_fit->SetParLimits(3, 0.5, 30.0);
    cb_fit->SetParLimits(4, 1.1, 100.0);
    cb_fit->SetLineColor(kRed);
    cb_fit->SetLineWidth(2);
    cb_fit->SetLineStyle(2);
    h1_Espec->Fit(cb_fit, "RQN");

    gaus_fit->Draw("same");
    cb_fit->Draw("same");

    auto *fit_leg = new TLegend(0.55, 0.55, 0.92, 0.9);
    fit_leg->SetFillStyle(0);
    fit_leg->SetBorderSize(0);
    fit_leg->AddEntry(h1_Espec, "Histogram", "l");

    TString gaus_label = Form("gaus: A=%.3f, #mu=%.3f, #sigma=%.3f",
                              gaus_fit->GetParameter(0),
                              gaus_fit->GetParameter(1),
                              gaus_fit->GetParameter(2));
    fit_leg->AddEntry(gaus_fit, gaus_label.Data(), "l");

    TString cb_label = Form("crystal ball: A=%.3f, #mu=%.3f, #sigma=%.3f, #alpha=%.3f, n=%.3f",
                            cb_fit->GetParameter(0),
                            cb_fit->GetParameter(1),
                            cb_fit->GetParameter(2),
                            cb_fit->GetParameter(3),
                            cb_fit->GetParameter(4));
    fit_leg->AddEntry(cb_fit, cb_label.Data(), "l");
    fit_leg->Draw();

    fit_compare->Update();

    if (file) file->Close();
    delete file;
}
