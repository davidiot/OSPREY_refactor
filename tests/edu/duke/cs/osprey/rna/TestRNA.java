package edu.duke.cs.osprey.rna;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.control.GMECFinder;
import org.junit.Test;

import static org.hamcrest.Matchers.is;
import static org.junit.Assert.assertThat;

/**
 * Unit tests pertaining to RNA functionality
 *
 * @author dz54
 */
public class TestRNA extends TestBase{

    @Test
    public void test1fxlFH() {

        String folder = "test/1fxlFH.junit/";

        String[] args = new String[]{"-c", folder + "KStar.cfg", "findGMEC",
                folder + "System.cfg", folder + "DEE.cfg"};

        ConfigFileParser cfp = new ConfigFileParser(args);//args 1, 3+ are configuration files
        cfp.loadData();

        GMECFinder gf = new GMECFinder();
        gf.init(cfp);
        ConfSearch.EnergiedConf gmec = gf.calcGMEC().get(0);
    }

}
