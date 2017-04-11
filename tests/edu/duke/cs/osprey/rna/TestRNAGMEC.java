package edu.duke.cs.osprey.rna;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.control.GMECFinder;
import org.junit.Test;


import static org.hamcrest.Matchers.is;
import static org.junit.Assert.assertThat;

/**
 * Unit tests pertaining to RNA GMEC functionality
 *
 * @author dzhou
 */
public class TestRNAGMEC extends TestBase{

    /**
     * Find a GMEC for an RNA double strand with some base flips
     */
    @Test
    public void test1cslHGMEC() {

        String folder = "test/1CSLH.junit/";

        String[] args = new String[]{"-c", folder + "KStar.cfg", "findGMEC",
                folder + "System.cfg", folder + "DEE.cfg"};

        ConfSearch.EnergiedConf gmec = checkGMEC(args);

        assertThat(gmec.getEnergy(), isRelatively(-40.963, 1e-3));
        // DZ: it looks like several assignments are pretty close in energy.
        // assertThat(gmec.getAssignments(), is(new int[] {30, 20, 48, 6}));
    }

    /**
     * Find a GMEC when both RNA and protein are involved
     * Note that the stack size might need to be increased.
     * I ran it with -Xss3000k -dzhou
     */
    @Test
    public void test1fxlFHGMEC() {

        String folder = "test/1fxlFH.junit/";

        String[] args = new String[]{"-c", folder + "KStar.cfg", "findGMEC",
                folder + "System.cfg", folder + "DEE.cfg"};

        ConfSearch.EnergiedConf gmec = checkGMEC(args);

        assertThat(gmec.getEnergy(), isRelatively(-120.530, 1e-3));
        // DZ: it looks like several assignments are pretty close in energy.
        // assertThat(gmec.getAssignments(), is(new int[] {9, 32, 0, 13, 13}));
    }

    /**
     * Tests DEEPer with RNA ring puckers
     * ran with -Xss10m
     */
    @Test
    public void test2OEUHDEEPer(){
        String folder = "test/2OEUH.deeper.junit/";
        String[] args = new String[]{"-c", folder + "KStar.cfg", "findGMEC",
                folder + "System.cfg", folder + "DEE.cfg"};

        ConfSearch.EnergiedConf gmec = checkGMEC(args);

        assertThat(gmec.getEnergy(), isRelatively(-62.173, 1e-3));
        assertThat(gmec.getAssignments(), is(new int[] {4, 6, 15}));
    }

    /**
     * Tests DEEPer with RNA ring puckers
     * ran with -Xss10m
     */
    @Test
    public void test3P49FHDEEPer(){
        String folder = "test/3P49FH.deeper.junit/";
        String[] args = new String[]{"-c", folder + "KStar.cfg", "findGMEC",
                folder + "System.cfg", folder + "DEE.cfg"};

        ConfSearch.EnergiedConf gmec = checkGMEC(args);

        assertThat(gmec.getEnergy(), isRelatively(-49.946, 1e-3));
        assertThat(gmec.getAssignments(), is(new int[] {5, 7, 36, 5}));
    }

    private ConfSearch.EnergiedConf checkGMEC(String[] args) {
        ConfigFileParser cfp = new ConfigFileParser(args);//args 1, 3+ are configuration files
        cfp.loadData();

        GMECFinder gf = new GMECFinder();
        gf.init(cfp);
        return gf.calcGMEC().get(0);
    }
}
