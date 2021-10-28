const int num_cores = 64; //number of cores in TIGRESS detector
Int_t num_known_sources = 4; //number of sources that can be used for calibration

//sorts an array into ascending order
//preserves correspondance between the main array and the 'carry' array
void insertion_sort(Double_t array[], Double_t carry[], int len) {
    int i;
    for (i = 1; i < len; i++) {
        Double_t val = array[i];
        Double_t carry_val = carry[i];
        for (int j = i - 1; j >= 0; j--) {
            if (array[j+1] < array[j]) {
                array[j+1] = array[j];
                carry[j+1] = carry[j];
                array[j] = val;
                carry[j] = carry_val;
            } else {
                break;
            }//else
        }//for
    }//for
}//insertion_sort

//loads the histograms from histogram root file
void load_histograms(const char histogram_filepath[], TH1F *hist[], Int_t source_count, Int_t num_cores) {
    TFile *infile = new TFile(histogram_filepath);
    TH2F * charge_matrix = (TH2F*)infile->Get(Form("source_%i", source_count));
    for (int i = 0; i < num_cores; i++) {
        //hist[i] = (TH1F*)infile->Get(Form("hist%i_%i",source_count,i));
        // create projection for each core
        hist[i] = (TH1F*)charge_matrix->ProjectionY(Form("hist%i_%i", source_count, i), i + 1, i + 1);
    }
}

//using gamma ray spectra and actual energy peaks, calculates a linear equation to calibrate the detector's outputs
void linear_calibration(TList *list, TH1F *hist[], Int_t num_peaks_used, Double_t energy[], Double_t energy_er[], Double_t gains[num_cores], Double_t offsets[num_cores]) {

    //arrays to store the centroids of the peaks
    Double_t centroids[num_cores][num_peaks_used];
    Double_t centroids_er[num_cores][num_peaks_used];

    for (int i = 0; i < num_cores; i++) {
        int max_bin = hist[i]->GetXaxis()->GetXmax();
        int min_bin = hist[i]->GetXaxis()->GetXmin();
        hist[i]->SetAxisRange(100, max_bin-100,"X");
        hist[i]->SetBinContent(1, 0);
        Double_t intr = hist[i]->Integral(min_bin, max_bin, "width");
        if (intr < 1000) {
            cout << i << " FAILED!" << endl;
            for (int j = 0; j < num_peaks_used; j++) {
                centroids[i][j] = -1;
            }//for
            continue;
        }//if

        //roughly locates the peaks in the spectrum
        TSpectrum *spec = new TSpectrum(2 * num_peaks_used);
        Int_t num_found = spec->Search(hist[i], 2, "", 0.5);

        cout << "Found " << num_found << " peaks in histogram." << endl;

        //if too many or too few peaks have been found,
        //something has gone wrong and we move on to the next core
        if (num_found != num_peaks_used) {
            for (int j = 0; j < num_peaks_used; j++) {
                centroids[i][j] = -1;
            }//for
            continue;
        }//if

        //we fit a gaussian distribution around each peak
        //using the TSpectrum's data as a starting point
        TF1 *fit[num_peaks_used];
        Double_t* x_pos_array = spec->GetPositionX();

        for (int j = 0; j < num_found; j++) {
            Double_t x_pos = x_pos_array[j];
            Int_t bin = hist[i]->GetXaxis()->FindBin(x_pos);
            Double_t y_pos = hist[i]->GetBinContent(bin);

            fit[j] = new TF1(Form("Fit %i-%i", i, j), "[0]*(exp(-((x-[1])^2/(2*[2]^2))))*(1+ROOT::Math::erf([5]*((x-[1]))/([2]*pow(2,0.5)))) + [3] + [4]*(0.5*(1-(ROOT::Math::erf(((x-[1])/([2]*2^(0.5)))))))", x_pos - 50, x_pos + 50);
            fit[j]->SetParameters(y_pos, x_pos, 1, 15, 1, -1);
            fit[j]->SetParLimits(0, 10, 1e6); //area
            fit[j]->SetParLimits(1, x_pos - 10, x_pos + 10); //centroid
            fit[j]->SetParLimits(2, 0.2, 15); //sigma
            fit[j]->SetParLimits(4, 0.1, 100); //magnitude of step in background noise
            fit[j]->SetParLimits(5, -10, -0.1); //background noise constant

            //fitting the equation and storing the calculated centroids
            hist[i]->Fit(fit[j], "RQ+");
            centroids[i][j] = fit[j]->GetParameter(1);
            centroids_er[i][j] = fit[j]->GetParError(1);
            cout << fit[j]->GetParameter(2) << endl;
            cout << "CENTROID PROCESSED: Graph " << i << " Guess: " << x_pos << " Actual: " << centroids[i][j] << endl;
        }//for

        //to make sure the centroids match up to the correct energies
        //we sort the centroids in ascending order
        //making sure the centroids_er keep the correspondance
        insertion_sort(centroids[i], centroids_er[i], num_peaks_used);
    }//for

    //we graph the centroids vs their corresponding energies
    //and fit a linear equation on the points
    for (int i = 0; i < num_cores; i++) {
        //if the histogram is empty, skip this core
        if (centroids[i][0] == -1) {
            gains[i] = -1;
            offsets[i] = -1;
            continue;
        }//if

        TGraphErrors *gr = new TGraphErrors(num_peaks_used, centroids[i], energy, centroids_er[i], energy_er);
        gr->Draw("AP");
        list->Add(gr); //adding the graph to be saved to root file later

        TF1 *coeffFit = new TF1("coeffFit", "[0] + [1]*x");
        gr->Fit(coeffFit, "Q+");
        //storing linear equation parameters
        gains[i] = coeffFit->GetParameter(1);
        offsets[i] = coeffFit->GetParameter(0);
    } //for
}//linear_calibration

//Main method to be executed by GRSISort
void linear_energy() {

    Double_t co60_ener[2] = {1173.240, 1332.508};
    Double_t co60_ener_e[2] = {0.003, 0.004};

    TList *list = new TList;

    TH1F *lin_hist[num_cores]; //histograms for each core

    //linear calibration equation parameters for each core
    Double_t lin_gains[num_cores];
    Double_t lin_offsets[num_cores];

    //getting calibration source
    string source[5];
    string rootfile[5];
    int iii = 0;
    int jjj = 0;
    string sourcefile = "sourcelist.dat";
    ifstream fp;
    fp.open(sourcefile.c_str());
    while (fp.good()) {
        fp >> source[iii] >> rootfile[iii];
        if(source[iii].compare("co60")==0) jjj=iii;
        iii++;
    }

    //loading in histograms from analysis root file
    load_histograms("charge_histograms.root", lin_hist, jjj, num_cores);

    //adding histograms to be saved later
    for (int i = 0; i < num_cores; i++) {
        if (lin_hist[i]->GetMaximum() != 0) {
            list->Add(lin_hist[i]);
        }//if
    }//for

    linear_calibration(list, lin_hist, 2, co60_ener, co60_ener_e, lin_gains, lin_offsets);
    //write out linear equations parameters
    ofstream coeff_fw("lin_energy_coeff.txt");

    for (int i = 0; i < num_cores; i++) {
        coeff_fw << lin_gains[i] << endl;
        coeff_fw << lin_offsets[i] << endl;
    }//for

    coeff_fw.close();

    //write out histograms and linear fit graphs for review
    TFile * outfile = new TFile("lin_cal_fits.root", "RECREATE");
    outfile->cd();
    list->Write();
    outfile->Close();
}//linear_energy_RF
