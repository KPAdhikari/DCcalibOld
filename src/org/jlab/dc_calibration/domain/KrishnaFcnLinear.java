/**
 *
 * @author KPAdhikari
 */
package org.jlab.dc_calibration.domain;

import static org.jlab.dc_calibration.domain.Constants.beta;
import static org.jlab.dc_calibration.domain.Constants.cos30;
import static org.jlab.dc_calibration.domain.Constants.nThBinsVz;
import static org.jlab.dc_calibration.domain.Constants.rad2deg;
import static org.jlab.dc_calibration.domain.Constants.thEdgeVzH;
import static org.jlab.dc_calibration.domain.Constants.thEdgeVzL;
import static org.jlab.dc_calibration.domain.Constants.wpdist;

import org.freehep.math.minuit.FCNBase;
import static org.jlab.dc_calibration.domain.Constants.nSectors;
import org.jlab.groot.data.GraphErrors;

public class KrishnaFcnLinear implements FCNBase {

    private int Sector;
    private int SL;
    private int ThBin;
    private GraphErrors[][][] profileX;
    private boolean isLinear = false;

    public KrishnaFcnLinear(int Sector, int SL, int ThBin, GraphErrors[][][] profileX, boolean isLinear) {
        this.Sector = Sector;
        this.SL = SL;
        this.ThBin = ThBin;
        this.profileX = profileX;
        this.isLinear = isLinear;
    }

    public double errorDef() {
        return 1;
    }

    public double valueOf(double[] par) {
        double delta = 0.;
        double chisq = 0.;
        double thetaDeg = 0.;
        double docaNorm = 0.;
        double measTime = 0.;
        double measTimeErr = 0.;
        double calcTime = 0.;

        for (int i = 0; i < profileX[Sector][SL][ThBin].getDataSize(0); i++) {
            docaNorm = profileX[Sector][SL][ThBin].getDataX(i);
            measTime = profileX[Sector][SL][ThBin].getDataY(i);
            measTimeErr = profileX[Sector][SL][ThBin].getDataEY(i);
            calcTime = calcTimeFunc(0, SL, docaNorm, par);

            if (measTimeErr == measTimeErr && measTimeErr > 0.0 && docaNorm < 0.9) {
                delta = (measTime - calcTime) / measTimeErr; // error weighted deviation
                chisq += delta * delta;
            }
        }

        return chisq;
    }

    protected double calcTimeFunc(int debug, int SL, double docaByDocaMax, double[] par) {
        double dMax = 2 * wpdist[SL];
        double x = docaByDocaMax * dMax;
        double v0Par = par[0];
        double calcTime = x / v0Par;
        return calcTime;
    }
}
