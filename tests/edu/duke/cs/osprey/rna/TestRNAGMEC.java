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
     * The four mutable residues are right between the flipped bases
     * This one will take a while to run
     */
    @Test
    public void test1cslHGMEC() {

        String folder = "test/1CSLH.junit/";

        String[] args = new String[]{"-c", folder + "KStar.cfg", "findGMEC",
                folder + "System.cfg", folder + "DEE.cfg"};

        ConfigFileParser cfp = new ConfigFileParser(args);//args 1, 3+ are configuration files
        cfp.loadData();

        GMECFinder gf = new GMECFinder();
        gf.init(cfp);
        ConfSearch.EnergiedConf gmec = gf.calcGMEC().get(0);

        assertThat(gmec.getEnergy(), isRelatively(-132.342, 1e-3));
        assertThat(gmec.getAssignments(), is(new int[] {4, 15, 35, 3}));
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

        ConfigFileParser cfp = new ConfigFileParser(args);//args 1, 3+ are configuration files
        cfp.loadData();

        GMECFinder gf = new GMECFinder();
        gf.init(cfp);
        ConfSearch.EnergiedConf gmec = gf.calcGMEC().get(0);

        assertThat(gmec.getEnergy(), isRelatively(-109.789, 1e-3));
        assertThat(gmec.getAssignments(), is(new int[] {3, 32, 0, 4, 13}));
    }

}
