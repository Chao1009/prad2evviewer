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
std::array<double, 3> fitGaus(TH1F *h, float expectPeak, TF1 **fitResult = nullptr)
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
    if (fitResult) *fitResult = new TF1(gaus);
    return {mean, sigma, chi2};
}

// Fit a peak near expectPeak with a Crystal Ball and return {mean, sigma, chi2/ndf}.
// Returns {0,0,0} if the histogram or fit is invalid.
std::array<double, 3> fitCrystalBall(TH1F *h, float expectPeak,
                                                  float alpha, float n,
                                                  TF1 **fitResult = nullptr)
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

    const double threshold = 0.05 * peakHeight;
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

    const double lo = h->GetBinCenter(std::max(1, leftBin));
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
    cb.SetParLimits(0, 0.0, 5. * peakHeight);
    cb.SetParLimits(1, lo, hi);
    cb.SetParLimits(2, 1e-6, std::max(hi - lo, 1e-3));
    cb.SetParLimits(3, 1.0, 5.0);
    cb.SetParLimits(4, 1.0, 20.0);

    int fitStatus = h->Fit(&cb, "RN");
    if (fitStatus != 0) {
        cb.SetParameters(peakHeight, peak0, sigma0, alpha, n);
        cb.FixParameter(3, alpha);
        cb.FixParameter(4, n);
        fitStatus = h->Fit(&cb, "RQN");
    }
    if (fitStatus != 0) return {0., 0., 0.};

    const double mean = cb.GetParameter(1);
    const double sigma = std::abs(cb.GetParameter(2));
    if (!std::isfinite(mean) || !std::isfinite(sigma) || sigma <= 0.)
        return {0., 0., 0.};

    const double chi2 = (cb.GetNDF() > 0) ? cb.GetChisquare() / cb.GetNDF() : 0.;
    if (fitResult) *fitResult = new TF1(cb);
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
    gaus->SetLineWidth(4);
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
        cb->SetLineWidth(4);
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

    TFile *file = TFile::Open("prad_025319_quick_check.root");
    TH1F *h1_Espec = nullptr;
    if (file && !file->IsZombie()) {
        h1_Espec = dynamic_cast<TH1F*>(file->Get("module_energy/h_W565"));
    }
    if (!h1_Espec) {
        h1_Espec = new TH1F("h1_Espec_demo", "Demo energy spectrum;Energy;Counts",
                            120, 0.0, 5000.0);
        for (int i = 0; i < 200000; ++i) {
            double x = gRandom->Gaus(2239.0, 80.0);
            if (x > 0.0) h1_Espec->Fill(x);
        }
    }

    h1_Espec->SetLineColor(kBlack);
    h1_Espec->SetLineWidth(4);
    h1_Espec->Draw("hist");

    const float expectedPeak = static_cast<float>(mu);
    const float cbAlpha = 1.5f;
    const float cbN = 5.0f;
    TF1 *gausFit = nullptr;
    TF1 *cbFit = nullptr;
    const auto gausResult = fitGaus(h1_Espec, expectedPeak, &gausFit);
    const auto cbResult = fitCrystalBall(h1_Espec, expectedPeak, cbAlpha, cbN, &cbFit);

    auto addFitLegendEntry = [](TLegend *legend, TF1 *fit, const char *label) {
        if (fit) legend->AddEntry(fit, label, "l");
    };

    if (gausFit) {
        gausFit->SetLineColor(kBlue + 1);
        gausFit->SetLineWidth(4);
        gausFit->Draw("same");
    }
    if (cbFit) {
        cbFit->SetLineColor(kRed + 1);
        cbFit->SetLineWidth(4);
        cbFit->SetLineStyle(2);
        cbFit->Draw("same");
    }

    auto *fitLegend = new TLegend(0.52, 0.62, 0.96, 0.92);
    fitLegend->SetFillStyle(0);
    fitLegend->SetBorderSize(0);
    fitLegend->AddEntry(h1_Espec, "Data", "l");
    if (gausFit) {
        const char *label = Form("Gaussian: A=%.3g, #mu=%.3f, #sigma=%.3f; #chi^{2}/ndf=%.3f",
                                 gausFit->GetParameter(0), gausFit->GetParameter(1),
                                 gausFit->GetParameter(2), gausResult[2]);
        addFitLegendEntry(fitLegend, gausFit, label);
    }
    if (cbFit) {
        const char *label = Form("Crystal Ball: A=%.3g, #mu=%.3f, #sigma=%.3f; #alpha=%.3f, n=%.3f; #chi^{2}/ndf=%.3f",
                                 cbFit->GetParameter(0), cbFit->GetParameter(1),
                                 cbFit->GetParameter(2), cbFit->GetParameter(3),
                                 cbFit->GetParameter(4), cbResult[2]);
        addFitLegendEntry(fitLegend, cbFit, label);
    }
    fitLegend->Draw();
    fit_compare->Update();
}
