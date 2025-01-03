/**
 * ChromHMM - automating chromatin state discovery and characterization 
 * Copyright (C) 2008-2012 Massachusetts Institute of Technology
 * 
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 **/
package edu.mit.compbio.ChromHMM;

import java.io.*;
import java.util.*;
import java.math.*;
import java.text.*;
import java.awt.*;
import java.awt.geom.*;
import java.awt.image.*;
import java.net.*;
import java.util.zip.*;
import org.tc33.jheatchart.HeatChart;

 
/**
 * The main class of ChromHMM implements command line parsing and core algorithms
 * The ChromHMM code was written by Jason Ernst 
 */
public class ChromHMM
{
  
    /**
     * True if should print out debug information
     */
    static boolean BVERBOSE = false;
    static double EPSILONEMISSIONS = Math.pow(10,-300);

    static double EPSILONSTATE = Math.pow(10,-300);

    static String SZSEGMENTEXTENSION = "_segments.bed";
    static String SZPOSTERIOREXTENSION = "_posterior.txt";
    static String SZSTATEBYLINEEXTENSION = "_statebyline.txt";

    /**
     * The default number of base pairs in a bin
     */
    static int DEFAULT_BINSIZEBASEPAIRS = 200;

    /**
     * When splitting the default number of bins per file
     */
    static int DEFAULT_NUMSPLITBINS = 5000;
 
    static int DEFAULT_OVERLAPENRICHMENT_NOFFSETLEFT = 0;
    static int DEFAULT_OVERLAPENRICHMENT_NOFFSETRIGHT = 1;
    static boolean DEFAULT_OVERLAPENRICHMENT_BCENTER = false;
    static boolean DEFAULT_OVERLAPENRICHMENT_BCOUNTMULTI = false;
    static boolean DEFAULT_OVERLAPENRICHMENT_BUSESIGNAL = false;
    static boolean DEFAULT_OVERLAPENRICHMENT_BBASERES = true;
    static boolean DEFAULT_OVERLAPENRICHMENT_BUNIFORMHEAT = false;

    static int DEFAULT_NEIGHBORHOOD_NUMLEFT = 10;
    static int DEFAULT_NEIGHBORHOOD_NUMRIGHT = 10;
    static boolean DEFAULT_NEIGHBORHOOD_BUSESTRAND = true; 
    static boolean DEFAULT_NEIGHBORHOOD_BUSESIGNAL= false;
    static int DEFAULT_NEIGHBORHOOD_NOFFSETANCHOR= 0;

    static String ANCHORFILEDIR = "ANCHORFILES";
    static String COORDDIR = "COORDS";
    static String CHROMSIZESDIR = "CHROMSIZES";

    static String SZNEIGHBORHOODEXTENSION = "_neighborhood";
    static String SZOVERLAPEXTENSION = "_overlap";
    static String SZBROWSERDENSEEXTENSION ="_dense";
    static String SZBROWSEREXPANDEDEXTENSION ="_expanded";

    /**
     * This is a parameter just to improve efficiency of the model learn does not impact results
     * in terms of when to try to exploit transition sparsity for improved effeciency
     */
    private static double SPARSECUTOFFRATIO = 0.7;

    /**
     * This is a parameter just to improve efficiency of the model learning does not impact results
     * in terms of when to try to exploit transition sparsity for improved effeciency
     */
    private static double SPARSECUTOFFLOOSERRATIO = 0.8;

    /**
     * Default Red value for heatmaps on 0 to 255 scale
     */
    static int DEFAULTCOLOR_R = 0;

    /**
     * Default Green value for heatmaps on 0 to 255 scale
     */
    static int DEFAULTCOLOR_G = 0;

    /**
     * Default Blue value for heatmaps on 0 to 255 scale
     */
    static int DEFAULTCOLOR_B = 255;

    static String[] ORDERSTRINGS = {"User","Emission","Transition","Fixed"};

    static char[] ORDERCHARS = {'U','E','T','F'};

    /**
     * Constant for user supplied state ordering
     */
    static int STATEORDER_USER = 0;

    /**
     * Constant for state ordering being based on the emission paraemters
     */
    static int STATEORDER_EMISSION = 1;

    /**
     * Constant for state ordering being based on the transition parameters
     */
    static int STATEORDER_TRANSITION = 2;

    /**
     *  Constant for state ordering being fixed
     */
    static int STATEORDER_FIXED = 3;


    /**
     * Constant for parameter initialization based on information measure
     */
    static int INITMETHOD_INFORMATION = 0;

    /**
     * Constant for parameter initialization being randomly selected
     */
    static int INITMETHOD_RANDOM = 1;

    /**
     * Constant for parameter initialization to be based on loading a model file     
     */
    static int INITMETHOD_LOAD = 2;

    /**
     * The log likelihood of the model
     */
    double dloglike;

    /**
     * Parameter used in the information based smoothing to smooth away from 0
     */
    double dinformationsmooth;


    /**
     * The directory where the model should be output
     */
    String szoutputdir = "";

    /**
     * The directory with the binarized input that should be used
     */
    String szinputdir = "";

    /**
     * Initial probability of each state
     */
    double[] probinit;

    /**
     * transitionprobs[ni][nj] contains the probability of transitioning from state i to state j
     */
    double[][] transitionprobs;

    /**
     * transitionprobs[ni] has the number of transitions from state i that have not been eliminated
     * because of low probability.
     */
    int[] transitionprobsnum;

    /**
     * transitionprobs[ni] from 0 to transitionprobsnum[ni]-1 has the indicies of the
     * states from state ni which have a non-eliminated transition
     */
    int[][] transitionprobsindex;

    /**
     * True if the transition between states has been eliminated, false otherwise.
     */
    boolean[][] elim;

    /**
     * Terminates if the estimated likelihood change after a full iteration is less than the value
     * If the value is negative this value is not used. The estimated likelihood is computed by adding
     * the likelihood of each sequence from its most recent pass through the genome.
     */
    double dconvergediff;


    /**
     * HashSet with all the file prefixes
     */
    HashSet hsprefix;

    /**
     * The maximum number of training iterations 
     */
    int nmaxiterations;

    /**
     * The number of non-eliminated states into state i
     */
    int[] transitionprobsnumCol;


    /**
     * trainsitionprobsindexCol[ni] from 0 to transitionprobsindexCol[ni]-1 has the indicies of the
     * states which have a non-eliminated transition
     */
    int[][] transitionprobsindexCol;


    /**
     * For each state, each feature we have an emission probability for each bucket
     */
    double[][][] emissionprobs;

    /**
     * The number of buckets for the emission parameter. 
     * The assumed data is binary but this is kept as a parameter for flexibility
     */
    private int numbuckets = 2;
 
    /**
     * The number of hidden states in the model 
     */
    int numstates;

    /**
     * traindingdataObservedIndex[ni][nj] indicates for the i^th sequence and j^th position
     * the combination which was observed.
     */
    int[][] traindataObservedIndex;

    /**
     * trainingdataObservedValues[ni][nj] indicates if for the i^th combination of marks the value
     * of the j^th mark
     */
    boolean[][] traindataObservedValues;

    /**
     * traindataNotMissing[ni][j] indicates if for the i^th combination of marks and missing values
     * the j^th mark is not missing
     */
     boolean[][] traindataNotMissing;

    /**
     * traindataObservedSeqFlags[ni][nj] contains whether on the i^th sequence the j^th combination of
     * marks were observed
     */
    boolean[][] traindataObservedSeqFlags;


    /**
     * The names of the data sets for which a joint model will be learned
     */ 
    String[] datasets;

    /**
     * The number of data sets that for which a joint model will be learned this is datasets.length
     */
    int numdatasets;


    /**
     * Stores the random number generator
     */
    Random theRandom;


    /**
     * Name of the file from which to base parameter initialization off of
     * If null then randomly initialize parameters
     */
    String szInitFile;


    /**
     * Determines how much weight is given to uniform vs. is the pre-loaded setting used to smooth around 0
     */
    double dloadsmoothtransition;

    /**
     * Determines how much weight is given to uniform vs. is the pre-loaded setting used to smooth around 0
     */
    double dloadsmoothemission;


    /**
     * Stores the cell associated with each sequence
     */
    String[] cellSeq;

    /**
     * Stores the chromosome associated with each sequence
     */
    String[] chromSeq;

    /**
     * File containing a list of 1-based state mapping
     */
    String szstateorderingfile;

    /**
     * File containing the desired order of columns listed sequentially
     */
    String szcolumnorderingfile;


    /**
     *  The set of chromosome sequence files from which to learn the model
     */
    String[] chromfiles;


    /**
     * Maps the internal state ordering to an actual ordering. This
     * is reflect in the emission, transition, and model files.
     */
    int[] stateordering;

    /**
     * Maps the column ordering to the ordering used in viewing the emission table
     * Note the columns are not reordered in the model file. They are the same as
     * the inital input files
     */
    int[] colordering;

    /**
     * Stores the selected parameter initialization method
     */
    int ninitmethod;

    /**
     * Stores the selected state ordering method
     */
    int nstateorder;

    /**
     * If true reorders columns in the emission matrix
     */
    boolean bordercols;


    /**
     * If true reorders rows of emission matrix
     */
    boolean borderrows;


    /**
     * File with list of input files to learn
     */
    String szinputfilelist;

    /**
     * If true prints to files the posterior distributions
     */    
    boolean bprintposterior;

    /**
     * If true prints a four column segmentation file
     */
    boolean bprintsegment;

    /**
     * If true prints the maximum state assignment one per line
     */
    boolean bprintstatebyline;


    /**
     * The number of base pairs in a bin
     */
    int nbinsize;

    /**
     * Transitions below 10^ntransitionpower will be eliminated for efficiency
     */
    int nzerotransitionpower;

    /**
     * Stores the Color to be used for the heatmap
     */
    Color theColor;


    /**
     * Character storing ordering method E-emission, T-transition, U-user
     */
    char chorder;

    /**
     * Descriptive string on the state ordering
     */
    String szorder;

    /**
     * An ID or name that will be included in some output files
     */
    String szoutfileID;

    /**
     * This parameter specifies the maximum number of seconds in learning and terminates
     * after a full iteration in which this is exceeded
     */
    int nmaxseconds;

    /**
     * Contains the maximum coordinate for each chromosome
     */
    String szchromlengthfile;

    /**
     * The header line of a loaded model
     */
    private String szLoadHeader;

    /**
     * Stores a mapping from states to labels
     */
    HashMap hmlabelExtend;

    ///////////////////////
    //code for confusion matrix
    /**
     * True iff EvalSubset should read input from posterior files
     */
    boolean breadposterior;

    /**
     * True iff EvalSubset should read input from standard segment files 
     */
    boolean breadsegment;

    /**
     * True iff EvalSubset should read input from segmentation file with one position per line
     */
    boolean breadstatebyline;

    /**
     * A bit string specifying for each mark whether each mark is included '1' or not included '0'
     */
    String szincludemarks;

    /**
     * True if the confusion matrix should be appended to the confusion output file
     */
    boolean bappend;

    /**
     * The directory containing the segmentations
     */
    String szsegmentdir;

    /**
     * The file to output confusion matrix 
     */
    String szconfusionfileprefix;

    /**
     * Variable which determines whether to do incremental or normal EM
     */
    boolean bnormalEM;

    /**
     * Maximum number of threads to try launching also constrained by what is available 
     * If less than one then set to maximum available
     */
    int nmaxprocessors;

    /**
     * stores length of each chromosome
     */
    int[] numtime;

    /**
     * whether to use low memory
     */
    boolean blowmem;

    /**
     * number of sequences to include in sampling
     */
    int numincludeseq;


    /**
     * Whether to print emission and transition parameter image files in addition to text files
     */
    boolean bprintimage;


    /**
     * true if emission parameters should be rescaled
     */
    boolean bscaleemissions = false;

    /**
     * true if should add pseudo-counts when learning models
     */
    boolean bpseudo = false;

    /**
     * true if should skip the last input from pseudo
     */
    boolean bnopseudolast = false;


    /**
     * true if output should be zipped files
     */
    boolean bgzip = false;


    /**
     * True if the binarized files aren't true chromosomes and instead are split portions
     */
    boolean bsplit = false;


    /**
     * The max number of bins per file when splitting 
     */
    int numsplitbins = ChromHMM.DEFAULT_NUMSPLITBINS;

    /**
     * BED file input that is remapped in Reorder
     */
    String szreorderinbedfile;

    /**
     * BED file output that was remapped in Reorder
     */
    String szreorderoutbedfile;


    /**
     * Indicate the beta's should be rescaled separately and also enforces that no beta value goes below STATEEPSILON or alpha if emission is greater than 0
     */
    boolean bscalebeta;

    /**
     * True if should also reorder the columns of the model file
     */
    boolean breordercolsmodel;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /**
     * Stores an integer index and array of boolean flags
     */
    static class ObservedRec
    {
       int nobserved;
       boolean[] flagA;

       ObservedRec(int nobserved, boolean[] flagA)
       {
	  this.nobserved = nobserved;
	  this.flagA = flagA;
       }
    }


    ///////////////////////////////////////////////////////////////////////////////////////////
    /**
     * Record stores an integer and random values, with reccompare can be used to
     * randomly sort a set of indicies
     */
    static class RecIntDouble
    {
	int nindex;
	double dval;

	RecIntDouble(int nindex, double dval)
	{
	    this.nindex = nindex;
	    this.dval = dval;
	}
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    /**
     * Sorts the records based on the value of dval
     */
    public static class RecIntDoubleCompare implements Comparator, Serializable
    {

	public int compare(Object o1, Object o2)
	{
	    RecIntDouble r1 = (RecIntDouble) o1;
	    RecIntDouble r2 = (RecIntDouble) o2;
	    if (r1.dval< r2.dval)
	    {
		return -1;
	    }
	    else if (r1.dval > r2.dval)
	    {
		return 1;
	    }
	    else
	    {
		return 0;
	    }
	}
    }

    ///////////////////////////////////////////////////////////////////////////////////////

    /**
     * Record stores an integer and random values, with reccompare can be used to
     * randomly sort a set of indicies
     */
    static class RecIntString
    {
	int nindex;
	String sz;

	RecIntString(int nindex, String sz)
	{
	    this.nindex = nindex;
	    this.sz = sz;
	}
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    /**
     * Sorts the records based on the value of dval
     */
    public static class RecIntStringCompare implements Comparator, Serializable
    {

	public int compare(Object o1, Object o2)
	{
	    RecIntString r1 = (RecIntString) o1;
	    RecIntString r2 = (RecIntString) o2;
	    return r1.sz.compareTo(r2.sz);
	}
    }


    ///////////////////////////////////////////////////////////////////////////////////////

    /**
     * Record stores the overall index, cell, chromosome values, and index within the chromosome which with RecIntStringSplitCompare can be used to
     * sort the files
     */
    static class RecIntStringSplit
    {
	int nindex; //overall index of file
	String szcell;
	String szchrom;
	int nsplitbinindex; //index within chromosome

	RecIntStringSplit(int nindex, String szcell, String szchrom, int nsplitbinindex)
	{
	    this.nindex = nindex;
	    this.szcell = szcell;
	    this.szchrom = szchrom;
	    this.nsplitbinindex = nsplitbinindex;
	}
    }


    ///////////////////////////////////////////////////////////////////////////////////////////
    /**
     * Sorts the records based on the value of szcell, then szchrom, then nsplitbinindex
     */
    public static class RecIntStringSplitCompare implements Comparator, Serializable
    {
	
	public int compare(Object o1, Object o2)
	{
	    //need to sort on chromosome name not file name
	    RecIntStringSplit r1 = (RecIntStringSplit) o1;
	    RecIntStringSplit r2 = (RecIntStringSplit) o2;

	    if (r1.szcell.equals(r2.szcell))
	    {
	       if (r1.szchrom.equals(r2.szchrom))
	       {
	 	  if (r1.nsplitbinindex < r2.nsplitbinindex)
		  {
		     return -1;
		  }
		  else if (r1.nsplitbinindex > r2.nsplitbinindex)
		  {
		     return 1;
		  }
		  else
		  {
		     return 0;
		  }
	       }
	       else
	       {
		   return (r1.szchrom+"_").compareTo(r2.szchrom+"_");
	       }
	    }
	    else
	    {
		return (r1.szcell+"_").compareTo(r2.szcell+"_");
	    }
	}
    }




    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * Constructor initializes the variable and loads the data used for learning the model
     */
     public ChromHMM(String szinputdir, String szoutputdir, String szinputfilelist,String szchromlengthfile, int numstates,int nseed, int ninitmethod,
                    String szInitFile, double dloadsmoothemission,double dloadsmoothtransition,double dinformationsmooth,
		     int nmaxiterations,double dcovergediff,int nmaxseconds,boolean bprintposterior,
		    boolean bprintsegment,boolean bprintstatebyline, int nbinsize,String szoutfileID,int nstateorder,boolean bordercols,int nzerotransitionpower,
		     Color theColor, boolean bnormalEM, int nmaxprocessors, boolean blowmem, 
                     int numincludeseq, boolean bprintimage, boolean bscaleemissions, 
                     boolean bpseudo, boolean bnopseudolast, boolean bgzip,boolean bsplit, boolean borderrows, boolean bscalebeta) throws IOException
    {
	this.szinputdir = szinputdir;
        this.szoutputdir = szoutputdir;
	this.szinputfilelist = szinputfilelist;
	this.szchromlengthfile = szchromlengthfile;
	this.numstates = numstates;
	this.ninitmethod = ninitmethod;
	this.szInitFile = szInitFile;
	this.dloadsmoothemission = dloadsmoothemission;
	this.dloadsmoothtransition = dloadsmoothtransition;
	this.dinformationsmooth = dinformationsmooth;
	this.nmaxiterations = nmaxiterations;
	this.nmaxseconds = nmaxseconds;
	this.dconvergediff = dcovergediff;
	this.bprintposterior = bprintposterior;
	this.bprintsegment = bprintsegment;
	this.bprintstatebyline = bprintstatebyline;
	this.nbinsize = nbinsize;
	this.szoutfileID = szoutfileID;
	this.nstateorder = nstateorder;
	this.chorder = ChromHMM.ORDERCHARS[nstateorder];
	this.szorder = ChromHMM.ORDERSTRINGS[nstateorder];
	this.bordercols = bordercols;
	this.borderrows = borderrows;
	this.nzerotransitionpower = nzerotransitionpower;
	this.theColor = theColor;
	this.bnormalEM = bnormalEM;
	this.nmaxprocessors = nmaxprocessors;
	this.numincludeseq = numincludeseq;
	this.blowmem = blowmem;
	this.bprintimage = bprintimage;
	this.bscaleemissions = bscaleemissions;
	this.bpseudo = bpseudo;
	this.bnopseudolast = bnopseudolast;
	this.bgzip = bgzip;
	this.bscalebeta = bscalebeta;
	//this.numsplitbins = numsplitbins;
	this.bsplit = bsplit;

        hmlabelExtend = new HashMap();
        theRandom = new Random(nseed);

	if (blowmem)
	{
	    loadDataFileStubs();
	}
	else
	{
	    loadData();
	}
	//loadData(); 
	//loadDataFileStubs();

	stateordering = new int[numstates];
	colordering = new int[numdatasets];
	for (int ni = 0; ni < stateordering.length; ni++)
	{
	    stateordering[ni] = ni;
	}

	for (int ni = 0; ni < colordering.length; ni++)
	{
	    colordering[ni] = ni;
	}
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////
    /**
     * Loads contents of szlabelmapping into a HashMap
     */
    private void makeLabelMapping(String szlabelmapping) throws IOException
    {
        if (szlabelmapping != null)
	{
	   BufferedReader bridlabel =  Util.getBufferedReader(szlabelmapping);
	   String szLine;

	    //Loading in a mapping from state ID to a label description

	    while ((szLine = bridlabel.readLine())!=null)
	    {
		StringTokenizer st = new StringTokenizer(szLine,"\t");
		String szID = st.nextToken().trim();
		String szLabelExtend = st.nextToken().trim();
		hmlabelExtend.put(szID,szLabelExtend);
	    }
	    bridlabel.close();
	}
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /**
     * Constructor for API use of ChromHMM that starts with a model loaded from a file
     */
    public ChromHMM(String szInitFile) throws IOException
    {

	this.szInitFile = szInitFile;
	loadModel();
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
     /**
      * Constructor involved in reordering the data
      */
    public ChromHMM(String szInitFile, String szoutputdir, String szstateorderingfile, String szcolumnorderingfile,
                    String szoutfileID,int nstateorder,boolean bordercols,Color theColor,
                    String szlabelmapping, boolean bprintimage,String szreorderinbedfile, String szreorderoutbedfile,
                    boolean breordercolsmodel) throws IOException
    {
	this.szcolumnorderingfile = szcolumnorderingfile;
	this.szInitFile = szInitFile;
	this.szoutputdir = szoutputdir;
	this.szstateorderingfile = szstateorderingfile;
	this.szoutfileID = szoutfileID;
	this.bprintimage = bprintimage;

	this.bordercols = bordercols;
	this.theColor = theColor;

	this.szreorderinbedfile = szreorderinbedfile;
        this.szreorderoutbedfile = szreorderoutbedfile;
	this.breordercolsmodel = breordercolsmodel;
        hmlabelExtend = new HashMap();
	makeLabelMapping(szlabelmapping);

	loadModel();

	//moved after loadModel in v1.12 since was getting overwritten if having user defined ordering 
	if (nstateorder != ChromHMM.STATEORDER_FIXED)
	{
	   this.nstateorder = nstateorder;
	}

	this.chorder = ChromHMM.ORDERCHARS[this.nstateorder];
	this.szorder = ChromHMM.ORDERSTRINGS[this.nstateorder];

	stateordering = new int[numstates];
	colordering = new int[numdatasets];
	for (int ni = 0; ni < stateordering.length; ni++)
	{
	    stateordering[ni] = ni;
	}

	for (int ni = 0; ni < colordering.length; ni++)
	{
	    colordering[ni] = ni;
	}
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /**
     * Constructor initializes the variable and loads the data used for making segmentation from a model
     */
    public ChromHMM(String szinputdir, String szinputfilelist, String szchromlengthfile, String szoutputdir, String szInitFile, String szoutfileID,
                    int nbinsize, boolean bprintposterior, boolean bprintsegment,boolean bprintstatebyline, 
                    boolean blowmem, boolean bscaleemissions, boolean bgzip, boolean bsplit, boolean bscalebeta) throws IOException
    {
	this.szinputdir = szinputdir;
	this.szinputfilelist = szinputfilelist;
	this.szchromlengthfile = szchromlengthfile;
	this.bprintposterior = bprintposterior;
	this.bprintsegment = bprintsegment;
	this.bprintstatebyline = bprintstatebyline;
	this.szoutfileID = szoutfileID;
        this.szoutputdir = szoutputdir;
	this.szInitFile = szInitFile;
	this.nbinsize = nbinsize;
	this.blowmem = blowmem;
	this.bscaleemissions = bscaleemissions;
	this.bgzip = bgzip;
	this.bsplit = bsplit;
	this.bscalebeta = bscalebeta;

        hmlabelExtend = new HashMap();

	if (blowmem)
	{
	   loadDataFileStubs();
	}
	else
        {
	   loadData();
	}

	loadModel();

	stateordering = new int[numstates];
	colordering = new int[numdatasets];
	for (int ni = 0; ni < stateordering.length; ni++)
	{
	    stateordering[ni] = ni;
	}

	for (int ni = 0; ni < colordering.length; ni++)
	{
	    colordering[ni] = ni;
	}
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /**
     * Constructor used for computing confusion results using a subset of marks for the EvalSubset command
     */
    public ChromHMM(String szinputdir, String szsegmentdir, String szinputfilelist, String szconfusionfileprefix,
                    String szInitFile, String szoutfileID,
                    int nbinsize, boolean breadposterior, boolean breadsegment,boolean breadstatebyline,
                    String szincludemarks, boolean bappend, Color theColor, boolean bprintimage, 
                    boolean blowmem, boolean bscaleemissions, boolean bscalebeta) throws IOException
    {
	this.bappend = bappend;
	this.szinputdir = szinputdir;
	this.szsegmentdir = szsegmentdir;
	this.szinputfilelist = szinputfilelist;
	//this.szchromlengthfile = szchromlengthfile; //removed since not doing anything
	this.breadposterior = breadposterior;
	this.breadsegment = breadsegment;
	this.breadstatebyline = breadstatebyline;
	this.szoutfileID = szoutfileID;
	this.szconfusionfileprefix = szconfusionfileprefix;
	this.szInitFile = szInitFile;
	this.nbinsize = nbinsize;
	this.szincludemarks = szincludemarks;
	this.theColor = theColor;
	this.bprintimage = bprintimage;
	this.blowmem = blowmem;
	this.bscaleemissions = bscaleemissions;
	this.bscalebeta = bscalebeta;
        hmlabelExtend = new HashMap();

	if (blowmem)
	{
	   loadDataFileStubs();
	}
	else
        {
	   loadData();
	}
	//loadData();
	loadModel();

	stateordering = new int[numstates];
	colordering = new int[numdatasets];
	for (int ni = 0; ni < stateordering.length; ni++)
	{
	    stateordering[ni] = ni;
	}

	for (int ni = 0; ni < colordering.length; ni++)
	{
	    colordering[ni] = ni;
	}
    }


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * Converts order abbreviation character to string
     */
    static String convertCharOrderToStringOrder(char ch)
    {
	for (int ni = 1; ni < ChromHMM.ORDERCHARS.length; ni++)
	{
	    if (ChromHMM.ORDERCHARS[ni] == ch)
	    {
		return ChromHMM.ORDERSTRINGS[ni];
	    }
	}
	return ChromHMM.ORDERSTRINGS[0];
    }



    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * Initializes the parameters based on the method determined by ninitmethod
     * then calls trainParameters
     */
    public void buildModel() throws IOException
    {
       if (ninitmethod ==ChromHMM.INITMETHOD_LOAD) 
       {
	   //loads parameters of the initial model
	   loadModelSmooth(dloadsmoothemission,dloadsmoothtransition);
       }
       else if (ninitmethod == ChromHMM.INITMETHOD_INFORMATION)
       {
	   if (blowmem)
	   {
	       informationInitializeNestedWithLoad();
	   }
	   else
	   {
	      informationInitializeNested();
	   }
       }
       else if (ninitmethod == ChromHMM.INITMETHOD_RANDOM)
       {
	   randomlyInitializeParams();
       }

       //trains the model
       //trainParameters();
       //trainParametersNormalEM();
       if (bnormalEM)
       {
	  if (blowmem)
	  {
             trainParametersParallelWithLoad();
	  }
	  else
	  {
	     trainParametersParallel();
	  }
       }
       else
       {
	   if (blowmem)
	   {
	      trainParametersWithLoad();
	   }
	   else
	   {
	      trainParameters();
	   }
       }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    private void reorderModel() throws IOException
    {	
       if (szstateorderingfile !=null)
       {
	   int[] stateorderingOLDtoNEW = new int[stateordering.length];
	   BufferedReader brstate =  Util.getBufferedReader(szstateorderingfile);
	   String szLine;
	   while ((szLine = brstate.readLine())!=null)
	   {
	       StringTokenizer st = new StringTokenizer(szLine,"\t ");
	       int nold = Integer.parseInt(st.nextToken())-1;
	       int nnew = Integer.parseInt(st.nextToken())-1;
	       stateordering[nnew] = nold;
	       stateorderingOLDtoNEW[nold] = nnew;
	   }
	   brstate.close();

	   if (szreorderinbedfile != null)
	   {
	       reorderBEDfile(stateorderingOLDtoNEW);
	   }
       }
       else
       {
          makeStateOrdering();
       }


       if (szcolumnorderingfile !=null)
       {
	   int ncol = 0;
	   BufferedReader brcol =  Util.getBufferedReader(szcolumnorderingfile);
	   String szLine;
	   HashMap hm = new HashMap();
	   while ((szLine = brcol.readLine())!=null)
	   {
	       hm.put(szLine,Integer.valueOf(ncol));
	       ncol++;
	   }
	   brcol.close();

	   for (int ni = 0;  ni < datasets.length; ni++)
	   {
	       Integer obj = ((Integer) hm.get(datasets[ni]));
	       if (obj == null)
	       {
		   throw new IllegalArgumentException(datasets[ni]+" not found!");
	       }
	       else
	       {
	          colordering[obj.intValue()] = ni;
	       }
	   }
       }
       else if (bordercols)
       {
          makeColOrdering();
       }
       //updates after each iteration the current status of the search
       printTransitionTable(-1);
       printEmissionTable(-1);

       if (bprintimage)
       {
          printEmissionImage(-1);
          printTransitionImage(-1);
       }
       printParametersToFile(-1);
    }


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    private void reorderBEDfile(int[] stateorderingOLDtoNEW) throws IOException
    {
	BufferedReader brin = Util.getBufferedReader(szreorderinbedfile);
	System.out.println("Writing to file "+szreorderoutbedfile);
	if (szreorderoutbedfile.toLowerCase(Locale.ENGLISH).endsWith(".gz"))
	{
	    GZIPOutputStream pwzip = new GZIPOutputStream(new FileOutputStream(szreorderoutbedfile));
	    //output zipfile
	    String szLine;
	    while ((szLine = brin.readLine())!=null)
	    {
		StringTokenizer st = new StringTokenizer(szLine, "\t");
		if (st.countTokens() < 4)
		{
		    //not full line just directly outputting except updating header info 
		    szLine = szLine.replaceAll("Emission ordered","User ordered");
		    szLine = szLine.replaceAll("Transition ordered","User ordered");
		    byte[] btformat = szLine.getBytes();
		    pwzip.write(btformat,0,btformat.length);
		    continue;
		}
		StringBuffer sb = new StringBuffer(st.nextToken().trim()+"\t"+st.nextToken().trim()+"\t"+st.nextToken().trim());
		String szinstate = st.nextToken().trim();
		String szoutstate;

		if (!Character.isDigit(szinstate.charAt(0)))
		{
		    szoutstate = (""+chorder) + (stateorderingOLDtoNEW[Integer.parseInt(szinstate.substring(1))-1]+1);
		}
		else
		{
                    szoutstate = ""+(stateorderingOLDtoNEW[Integer.parseInt(szinstate)-1]+1);
		}
		sb.append("\t"+szoutstate);

		while (st.hasMoreTokens())
		{
		    sb.append("\t"+st.nextToken().trim());
		}
		sb.append("\n");

		byte[] btformat = sb.toString().getBytes();
		pwzip.write(btformat,0,btformat.length);
	    }
	    pwzip.finish();
	    pwzip.close();
	    brin.close();
	}
	else
        {
	    PrintWriter pw = new PrintWriter(new FileWriter(szreorderoutbedfile));
	    //output zipfile
	    String szLine;
	    while ((szLine = brin.readLine())!=null)
	    {
		StringTokenizer st = new StringTokenizer(szLine, "\t");
		if (st.countTokens() < 4)
		{
		    //not full line just directly outputting except updating header info
		    szLine = szLine.replaceAll("Emission ordered","User ordered");
                    szLine = szLine.replaceAll("Transition ordered","User ordered");
		    pw.println(szLine);
		    continue;
		}
		pw.print(st.nextToken().trim()+"\t"+st.nextToken().trim()+"\t"+st.nextToken().trim());
		String szinstate = st.nextToken().trim();
		String szoutstate;
		if (!Character.isDigit(szinstate.charAt(0)))
		{
		    szoutstate = (""+chorder) + (stateorderingOLDtoNEW[Integer.parseInt(szinstate.substring(1))-1]+1);
		}
		else
		{
                    szoutstate = ""+ (stateorderingOLDtoNEW[Integer.parseInt(szinstate)-1]+1);
		}
		pw.print("\t"+szoutstate);

		while (st.hasMoreTokens())
		{
		    pw.print("\t"+st.nextToken().trim());
		}
		pw.println();
	    }
	    brin.close();
	    pw.close();
	}
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /**
     * If nstateorder is ChromHMM.STATEORDER_EMISSION orders emission parameters based on emission parameters and stores in stateordering
     * If nstateorder is ChromHMM.STATEORDER_TRANSITION orders emission parameters based on transition parameters and stores in stateordering
     */
    private void makeStateOrdering()
    {
	if (nstateorder ==  ChromHMM.STATEORDER_EMISSION)
	{
	   double[][] emissionprobspos = new double[numstates][numdatasets];
	   //just gets out the emissionprobability in index 1 of the third index
	   for (int ni = 0; ni < numstates; ni++)
	   {
	      double[] emissionprobspos_ni = emissionprobspos[ni];
	      double[][] emissionprobs_ni = emissionprobs[ni];
	      for (int nj = 0; nj < numdatasets; nj++)
	      {
	         emissionprobspos_ni[nj] = emissionprobs_ni[nj][1];
	      }
	   }

	   makeOrderingCorrelation(emissionprobspos, stateordering);
	}
	else if (nstateorder == ChromHMM.STATEORDER_TRANSITION)
	{
	    makeOrderingTransition(stateordering);
	}
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /**
     * Stores in colordering an ordering of the columns of the emission parameter matrix based on parameter correlaton
     */
    public void makeColOrdering()
    {

	double[][] emissionprobspostranspose = new double[numdatasets][numstates];
	for (int ni = 0; ni < numdatasets; ni++)
	{
	    double[] emissionprobspostranspose_ni = emissionprobspostranspose[ni];
	    for (int nj = 0; nj < numstates; nj++)
	    {
		emissionprobspostranspose_ni[nj] = emissionprobs[nj][ni][1];
	    }
	}
	makeOrderingCorrelation(emissionprobspostranspose,colordering);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /**
     * Returns an ordering of states using an approximation algorithm to minimize the total distance between neighboring states 
     * where distance between neighboring states (i,j) is defined as 2- (t_i,j +t_j,i)
     * The ordering is the best greedy ordering considering each state as the initial state
     */
    private void makeOrderingTransition(int[] ordering)
    {
	boolean[] assignedrow = new boolean[ordering.length];
        int[] temproworder = new int[ordering.length];
        double dmintotsum = Double.MAX_VALUE;

        for (int ninitrow = 0; ninitrow < ordering.length; ninitrow++)
	{
	    //tries ordering starting from each possible initial location
	    temproworder[0] = ninitrow;

	    for (int ni = 0; ni < assignedrow.length; ni++)
	    {
		assignedrow[ni] = false;
	    }

	    assignedrow[ninitrow] = true;

	    double dtotsum = 0;
	    int nminrow;
	    int nprevminrow=ninitrow;

	    for (int ncurrow = 1; ncurrow < ordering.length; ncurrow++)
	    {
		//finding the next state for the ordering
		double dmindist = Double.MAX_VALUE;
		nminrow = 0;
		//considering all not already assigned states
		for (int nrow = 0; nrow < ordering.length; nrow++)
		{
		   if (!assignedrow[nrow])
	           {
		       double ddist = 2-(transitionprobs[nrow][nprevminrow]+transitionprobs[nprevminrow][nrow]);
		       //checks if distance is less than minimum distance found so far, if so uses it
		       if (ddist < dmindist)
		       {
			   dmindist= ddist;
			   nminrow = nrow;
		       }
		   }
		}

		dtotsum += dmindist;  //increment total sum
		temproworder[ncurrow] = nminrow;
		assignedrow[nminrow] = true;
		nprevminrow = nminrow; //this is now the best previous row
	    }

	    if (dtotsum < dmintotsum)
	    {
		//best one found so far updating totalsum and storing roworder
		dmintotsum = dtotsum;
	        for (int ni = 0; ni < ordering.length; ni++)
	      	{
	       	    ordering[ni] = temproworder[ni];
		}
	    }
	}
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /**
     * Returns an ordering of states using an approximation algorithm to minimize the total distance between neighboring states 
     * where distance between neighboring states (i,j) is defined as the sqrt of 1-minus the correlation coefficient
     * The ordering is the best greedy ordering considering each state as the initial state
     */
    private void makeOrderingCorrelation(double[][] data, int[] ordering)
    {
	boolean[] assignedrow = new boolean[ordering.length];
        int[] temproworder = new int[ordering.length];
        double dmintotsum = Double.MAX_VALUE;

	double[][] correlationdistance = new double[ordering.length][ordering.length];
	for (int ni = 0; ni < ordering.length; ni++)
	{
	    for (int nj = 0; nj < ordering.length; nj++)
	    {
		correlationdistance[ni][nj] = Math.sqrt(1-Util.correlation(data[ni],data[nj]));
	    }
	}

        for (int ninitrow = 0; ninitrow < ordering.length; ninitrow++)
	{
	    temproworder[0] = ninitrow;

	    for (int ni = 0; ni < assignedrow.length; ni++)
	    {
		assignedrow[ni] = false;
	    }
	    assignedrow[ninitrow] = true;

	    double dtotsum = 0;
	    int nminrow;
	    int nprevminrow = ninitrow;

	    for (int ncurrow = 1; ncurrow < ordering.length; ncurrow++)
	    {
		//finding the next state for the ordering
		double dmindist = Double.MAX_VALUE;
		nminrow = 0;
		//considering all not already assigned states
		for (int nrow = 0; nrow < ordering.length; nrow++)
		{
		   if (!assignedrow[nrow])
	           {
		       double ddist = correlationdistance[nrow][nprevminrow];
		       //checks if distance is less than minimum distance found so far, if so uses it
		       if (ddist < dmindist)
		       {
			   dmindist= ddist;
			   nminrow = nrow;
		       }
		   }
		}

		dtotsum += dmindist; //increment total sum
		temproworder[ncurrow] = nminrow;
		assignedrow[nminrow] = true;
		nprevminrow = nminrow; //this is now the best previous row
	    }
	    if (dtotsum < dmintotsum)
	    {
		//best one found so far updating totalsum and storing roworder
		dmintotsum = dtotsum;
	        for (int ni = 0; ni < ordering.length; ni++)
	      	{
	       	    ordering[ni] = temproworder[ni];
		}
	    }
	}
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    /**
     * Outputs the transition parameters as a '.png' and 'svg' heatmap in the szoutputdir directory
     * the file name startswith 'transitions_numstates' and if szoutfileID is non-empty that is
     * included as well.
     */
    public void printTransitionImage(int niteration) throws IOException
    {
	//stores in sorteddata the contents of transitionprobs with the states 
	//reordered based on stateordering
	double[][] sorteddata = new double[numstates][numstates];
	String[] statelabels =  new String[numstates];

        for (int ni = 0; ni < numstates; ni++)
	{
	   double[] sorteddata_ni = sorteddata[ni];
	   double[] transitionprobs_stateordering_ni = transitionprobs[stateordering[ni]];
	   for (int nj =0; nj < numstates; nj++)
	   {
	      sorteddata_ni[nj] = transitionprobs_stateordering_ni[stateordering[nj]];
           }
        }

	for (int ni = 0; ni < numstates; ni++)
	{
	    statelabels[ni] = ""+(ni+1);
	    String szsuffix;
	    //changing in v1.12 to be state labels for new state to match char ordering
	    //if ((szsuffix = (String) hmlabelExtend.get(""+chorder+(stateordering[ni]+1)))!=null)
	    if ((szsuffix = (String) hmlabelExtend.get(""+chorder+(ni+1)))!=null)
	    {
	       statelabels[ni]+= "_"+szsuffix;
	    }
	}

        HeatChart map = new HeatChart(sorteddata);

        map.setTitle("Transition Parameters");
        map.setXAxisLabel("State To ("+szorder+" order)");
        map.setYAxisLabel("State From ("+szorder+" order)");
	map.setAxisValuesFont(new Font("SansSerif",0,20));
	map.setAxisLabelsFont(new Font("SansSerif",0,22));
	map.setTitleFont(new Font("SansSerif",0,24));
	if (sorteddata.length <=5)
	{
	   map.setChartMargin(125);
	}
	else
	{
	   map.setChartMargin(50);
	}
        map.setXValues(statelabels);  //sets the state labels to be the state IDs
        map.setYValues(statelabels);
        map.setHighValueColour(theColor);
 
	String szfileprefix; 
	if (szoutfileID.equals(""))
	{
	   szfileprefix = szoutputdir+"/transitions_"+numstates;
	}
	else
	{
	   szfileprefix = szoutputdir+"/transitions_"+numstates+"_"+szoutfileID;
	}

        map.saveToFile(new File(szfileprefix+".png"));
	Util.printImageToSVG(map, szfileprefix+".svg");

	if (niteration <= 1)
	{
	   System.out.println("Writing to file "+szfileprefix+".png");
	   System.out.println("Writing to file "+szfileprefix+".svg");
	}
    }

    ///////////////////////////////////////////////////////////////////////////////////////////

    /**
     * Outputs the emission parameters as a '.png' and '.svg' heatmap in the szoutputdir directory
     * the file name startswith 'emissions_numstates' and if szoutfileID is non-empty that is
     * included as well.
     */
    public void printEmissionImage(int niteration) throws IOException
    {
	double[][] sorteddata =new double[numstates][numdatasets];
	String[] sortedcollabels = new String[numdatasets];
	String[] rowlabels = new String[numstates];

	//copying in the emission parameters for the positive on bucket
        for (int ni = 0; ni < sorteddata.length; ni++)
	{
	   for (int nj =0; nj < sorteddata[ni].length; nj++)
	   {
	      	sorteddata[ni][nj] = emissionprobs[stateordering[ni]][colordering[nj]][1];
           }
        }

        for (int ni = 0; ni < numdatasets; ni++)
        {
	    //column labels also need to be reordered
           sortedcollabels[ni] = datasets[colordering[ni]];
	}

	for (int ni = 0; ni < numstates; ni++)
	{
	    rowlabels[ni] = ""+(ni+1);
	    String szsuffix;
	    //changing in v1.12 to be state labels for new state to match char ordering
	    //if ((szsuffix = (String) hmlabelExtend.get(""+chorder+(stateordering[ni]+1)))!=null)
	    if ((szsuffix = (String) hmlabelExtend.get(""+chorder+(ni+1)))!=null)
	    {
	       rowlabels[ni]+= "_"+szsuffix;
	    }
	}

        HeatChart map = new HeatChart(sorteddata);

        map.setTitle("Emission Parameters");
        map.setXAxisLabel("Mark");
	map.setAxisValuesFont(new Font("SansSerif",0,20));
	map.setAxisLabelsFont(new Font("SansSerif",0,22));
	map.setTitleFont(new Font("SansSerif",0,24));
        map.setYAxisLabel("State ("+szorder+" order)");
	if (sorteddata.length <=5)
	{
	   map.setChartMargin(125);
	}
	else
	{
	   map.setChartMargin(50);
	}
        map.setXValues(sortedcollabels);
        map.setYValues(rowlabels);
        map.setHighValueColour(theColor);

	String szfileprefix;
	if (szoutfileID.equals(""))
	{
	   szfileprefix = szoutputdir+"/emissions_"+numstates;
	}
	else
	{
	   szfileprefix = szoutputdir+"/emissions_"+numstates+"_"+szoutfileID;
	}

	Util.printImageToSVG(map, szfileprefix+".svg");
	map.saveToFile(new File(szfileprefix+".png"));

	if (niteration <= 1)
	{
	   System.out.println("Writing to file "+szfileprefix+".svg");
	   System.out.println("Writing to file "+szfileprefix+".png");
	}
    }

    ///////////////////////////////////////////////////////////////////////////////////////////

    /**
     * Prints the contents of the mark present emission parameters out to a text file with the 
     * states and marks ordered based on stateordering and colordering
     */
    public void printEmissionTable(int niteration) throws IOException
    {
        PrintWriter pw;

	String szfile;
	if (szoutfileID.equals(""))
	{
	    szfile = szoutputdir+"/emissions_"+numstates+".txt";
           pw = new PrintWriter(szfile);
	}
	else
	{
	   szfile = szoutputdir+"/emissions_"+numstates+"_"+szoutfileID+".txt";
           pw = new PrintWriter(szfile);
	}

	if (niteration <= 1)
	{
           System.out.println("Writing to file "+szfile);
	}

	pw.print("State ("+szorder+" order)");
	for (int ni = 0; ni < datasets.length; ni++)
	{
	    pw.print("\t"+datasets[colordering[ni]]);
	}
	pw.println();


	for (int ni = 0; ni < numstates; ni++)
	{
	   pw.print((ni+1));
           String szsuffix;
	    //changing in v1.12 to be state labels for new state to match char ordering
           //if ((szsuffix = (String) hmlabelExtend.get(""+chorder+(stateordering[ni]+1)))!=null)
           if ((szsuffix = (String) hmlabelExtend.get(""+chorder+(ni+1)))!=null)
	   {
	       pw.print("_"+szsuffix);
	   }

           for (int nj = 0; nj < emissionprobs[ni].length; nj++)
	   {
	      pw.print("\t"+emissionprobs[stateordering[ni]][colordering[nj]][1]);
	   }
	   pw.println();
	}
	pw.close();
    }

    ///////////////////////////////////////////////////////////////////////////////////////////

    /**
     *  Prints the contents of the transition parameters out to a text file with the states ordered
     */
    public void printTransitionTable(int niteration) throws IOException
    {
	PrintWriter pw;
	String szfile;
	if (szoutfileID.equals(""))
	{
	    szfile = szoutputdir+"/transitions_"+numstates+".txt";
	    pw = new PrintWriter(szfile);
	}
	else
	{
	    szfile = szoutputdir+"/transitions_"+numstates+"_"+szoutfileID+".txt";
	    pw = new PrintWriter(szfile);
	}

	if (niteration <= 1)
	{
	   System.out.println("Writing to file "+szfile);
	}

	pw.print("State (from\\to) ("+szorder+" order)");
	for (int ni = 0; ni < numstates; ni++)
	{
	    pw.print("\t"+(ni+1));
	    String szsuffix;
	    //changing in v1.12 to be state labels for new state to match char ordering
	    //if ((szsuffix = (String) hmlabelExtend.get(""+chorder+(stateordering[ni]+1)))!=null)
	    if ((szsuffix = (String) hmlabelExtend.get(""+chorder+(ni+1)))!=null)
	    {
	       pw.print("_"+szsuffix);
	    }
	}
	pw.println();


	for (int ni = 0; ni < numstates; ni++)
	{
	    pw.print(""+(ni+1));
	    String szsuffix;
	    //changing in v1.12 to be state labels for new state to match char ordering
	    //if ((szsuffix = (String) hmlabelExtend.get(""+chorder+(stateordering[ni]+1)))!=null)
	    if ((szsuffix = (String) hmlabelExtend.get(""+chorder+(ni+1)))!=null)
	    {
	       pw.print("_"+szsuffix);
	    }
	   double[] transitionprobs_stateordering_ni = transitionprobs[stateordering[ni]];
           for (int nj = 0; nj < numstates; nj++)
	   {
	      pw.print("\t"+transitionprobs_stateordering_ni[stateordering[nj]]);
	   }
	   pw.println();
	}
	pw.close();
    }

    ///////////////////////////////////////////////////////////////////////////////////////////

    /**
     * Prints to a file the contents of the model. States are printed based on state-ordering
     * but columns are not reordered and are consistent as in the original data file
     * szoutputdir+"/model_"+numstates+"_"+szoutfileID+".txt"
     * if -1 then uses model header
     */
    private void printParametersToFile(int niteration) throws IOException
    {
	PrintWriter pw;

	String szfile;
	if (szoutfileID.equals(""))
	{
	  szfile = szoutputdir+"/model_"+numstates+".txt";
          pw = new PrintWriter(szfile);
	}
	else
	{
	   szfile = szoutputdir+"/model_"+numstates+"_"+szoutfileID+".txt";
           pw = new PrintWriter(szfile);
	}

	if (niteration <= 1)
	{
	   System.out.println("Writing to file "+szfile);
	}
	
	if (niteration == -1)
	{
	    StringTokenizer st = new StringTokenizer(szLoadHeader);
	    pw.print(st.nextToken()+"\t"+st.nextToken()+"\t"+chorder);
	    st.nextToken(); //might be changing the order type
	    if (st.hasMoreTokens())
	    {
		//added this check in version 1.11 since elim models won't have these two tokens for loglikelihood and iteration
	       pw.println("\t"+st.nextToken()+"\t"+st.nextToken());
	    }
	}
	else
	{
	   pw.println(numstates+"\t"+numdatasets+"\t"+chorder+"\t"+dloglike+"\t"+niteration);
	}

	for (int ni = 0; ni < numstates; ni++)
	{
	    pw.println("probinit\t"+(ni+1)+"\t"+probinit[stateordering[ni]]);
	}

	for (int ni = 0; ni < transitionprobs.length; ni++)
        {
	   double[] transitionprobs_stateordering_ni = transitionprobs[stateordering[ni]];
           for (int nj = 0; nj < transitionprobs_stateordering_ni.length; nj++)
	   {
	       pw.println("transitionprobs\t"+(ni+1)+"\t"+(nj+1)+"\t"+transitionprobs_stateordering_ni[stateordering[nj]]);
	   }
	}

	for (int ni = 0; ni < numstates; ni++)
	{
	   double[][] emissionprobs_stateordering_ni = emissionprobs[stateordering[ni]];
           for (int nj = 0; nj < emissionprobs_stateordering_ni.length; nj++)
	   {
	      double[] emissionprobs_stateordering_ni_nj;
	      String szdataset;
	      if (breordercolsmodel)
	      {
	         emissionprobs_stateordering_ni_nj = emissionprobs_stateordering_ni[colordering[nj]];
	         szdataset = datasets[colordering[nj]];
	      }
	      else
	      {
	         emissionprobs_stateordering_ni_nj = emissionprobs_stateordering_ni[nj];
	         szdataset = datasets[nj];
	      }

	      for (int nk = 0; nk < emissionprobs_stateordering_ni_nj.length; nk++)
	      {
	          pw.println("emissionprobs\t"+(ni+1)+"\t"+nj+"\t"+szdataset+"\t"+nk+"\t"+emissionprobs_stateordering_ni_nj[nk]);	  
	      }
	   }
	}
	pw.close();
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////




    public void informationInitializeNestedWithLoad() throws IOException
    {
	//inital probability vector
        probinit = new double[numstates];

	//creates the emission probability matrix
	emissionprobs = new double[numstates][numdatasets][numbuckets];

	//set to true if a transition has been eliminated
	elim = new boolean[numstates][numstates];

	//initalize the transition matrix
	transitionprobs = new double[numstates][numstates];

	//initialize index of the next non-zero transition
	transitionprobsindex = new int[numstates][numstates];

	//initalize number of non-zero transitions
	transitionprobsnum = new int[numstates];

	//initalize column-wise index of non-zero transitions
        transitionprobsindexCol = new int[numstates][numstates];
 
	//number of non-zero column transitions
        transitionprobsnumCol = new int[numstates];

	int nmaxtime = 0;
	for (int nseq = 0; nseq < numtime.length; nseq++)
	{
	   //numtime[nseq] = traindataObservedIndex[nseq].length;
	   if (numtime[nseq] > nmaxtime)
	   {
	      nmaxtime = numtime[nseq];
	   }
	}



	int[][] traindataObservedIndexPair = new int[chromfiles.length][];

	//stores the indicies of the data
	int[] traindataObservedIndex = new int[nmaxtime];

	//saving the mapping of signatures and chromsome observed on
	//stores whether there is a present call at each location
	boolean[][] traindataObservedValues = new boolean[nmaxtime][numdatasets];

	//stores whether the mark is not considered missing
	boolean[][] traindataNotMissing = new boolean[nmaxtime][numdatasets]; //usually nobserved
	

	ArrayList alobservedpairflags = new ArrayList(); //is an index from the element combination to the associated flags

	HashMap hmObserved= new HashMap(); //is an index from the associated flags to the element index

	int nobserved= 0; //index on the unique observation combination we are observing

        int[] samples= null;

        int ncurrnumincludeseq;
        if (numincludeseq >= 1)
        {
	   ncurrnumincludeseq = Math.min(numincludeseq, chromfiles.length);
	   samples = new int[ncurrnumincludeseq];
        }
        else
        {
	   ncurrnumincludeseq = chromfiles.length;
        }

        boolean[] bincludeseq = new boolean[chromfiles.length];
        for (int nk = 0; nk < bincludeseq.length; nk++)
        {
	   bincludeseq[nk] = true;
        }

        if (numincludeseq >=1)
        {
	   for (int ni = 0; ni < samples.length; ni++)
	   {
              samples[ni] = ni;
	   }

	   //Random theRandom2 = new Random(412);//433

	   //sampling without replacement data locations wanted
	   for (int ni = ncurrnumincludeseq; ni < chromfiles.length; ni++)
           {
	      if (theRandom.nextDouble() < ncurrnumincludeseq/((double) ni +1))
	      {
      	         samples[theRandom.nextInt(ncurrnumincludeseq)] = ni;
      	      }
	   }

           for (int nk = 0; nk < bincludeseq.length; nk++)
           {
              bincludeseq[nk] = false;
	   }

           for (int nk = 0; nk < samples.length; nk++)
	   {
	      bincludeseq[samples[nk]] = true;
	      //System.out.println("including\t"+samples[nk]+"\t"+bincludeseq[samples[nk]]);
           }
	}

	//generates all vectors of combinations of consecutive of '1' calls
	for (int nseq = 0; nseq < chromfiles.length; nseq++)
	{
	    if (bincludeseq[nseq])
	    {
	       //going through each sequence

	       //int[] traindataObservedIndex_nseq = traindataObservedIndex[nseq];
	       //int traindataObservedIndex_nseq_m1 = traindataObservedIndex_nseq.length -1;
	       int numtime_nseq_m1 = numtime[nseq] - 1;

	       HashMap hmObservedLoad = new HashMap(); //maps an observation string to an index and set of flags
 
    	       int nobservedload = 0;

	       if (ChromHMM.BVERBOSE)
	       {
                  System.out.println("reading\t"+szinputdir+" "+chromfiles[nseq]);
	       }
	      
               BufferedReader br = Util.getBufferedReader(szinputdir+"/"+chromfiles[nseq]);
               String szLine = br.readLine(); //first line tells cell type and chromosome
	       br.readLine();//flush mark header
	       ArrayList aldata = new ArrayList();
	       while ((szLine = br.readLine())!=null)
	       {
	          StringTokenizer st = new StringTokenizer(szLine,"\t ");
      	          StringBuffer sb = new StringBuffer();
		
	          for (int ncol = 0; ncol < numdatasets; ncol++)
      	          {
	             if (!st.hasMoreTokens())
	             {
	                throw new IllegalArgumentException("Found line without "+numdatasets+" values in file "+chromfiles[nseq]);
		     }

		     String sztoken = st.nextToken();//.trim();
		    
		     if (sztoken.equals("0"))
	             {
	                sb.append("0");
		     }
	             else if (sztoken.equals("1"))
	             {
	                sb.append("1");
		     }
	             else if (sztoken.equals("2"))
		     {
	                //this means missing
	                sb.append("2");
		     }
	             else
	             {
                        throw new IllegalArgumentException("Unrecognized value "+sztoken+" found in "+szinputdir+"/"+chromfiles[nseq]);
		     }
		  }
	          aldata.add(sb.toString());
	       }
	       br.close();
	 
	       int nsize = aldata.size();

	      //traindataObservedIndex[nfile] = new int[nsize];
	      //int[] traindataObservedIndex_nfile = traindataObservedIndex[nfile];

	       for (int nrow = 0; nrow < nsize; nrow++)
	       {
	          BigInteger theBigInteger = new BigInteger((String) aldata.get(nrow),3);
		  Integer theObservedInt  = (Integer) hmObservedLoad.get(theBigInteger);
	          //boolean[] flagA;

	          if (theObservedInt == null)
       	          {
		    //System.out.println(szmappingbyte.length());
		    //storing a mapping from observed byte string to an integer index in alFlags and alObserved
		      hmObservedLoad.put(theBigInteger, Integer.valueOf(nobservedload));

		     //saving this observed index
		     traindataObservedIndex[nrow] = nobservedload;

		     //increments the number of observed combinations of marks
		     nobservedload++;
		  }
		  else
		  {
		     //storing the index of the flags associated with this row 
		     traindataObservedIndex[nrow] = ((Integer) theObservedInt).intValue();
		  }
	       }	 
	    
	       Iterator hmObservedIterator = hmObservedLoad.entrySet().iterator();
	       while (hmObservedIterator.hasNext())
	       {
	          Map.Entry pairs = (Map.Entry) hmObservedIterator.next();
	          BigInteger theBigInteger = (BigInteger) pairs.getKey();
	          String szmapping = theBigInteger.toString(3);  //getting back the mapping string

	          //ObservedRec theObservedRec = (ObservedRec) pairs.getValue();
	          int ncurrindex = ((Integer) pairs.getValue()).intValue();// theObservedRec.nobserved;//this is an index on which obervation combination it is

	          boolean[] traindataObservedValues_ncurrindex = traindataObservedValues[ncurrindex];
	          boolean[] traindataNotMissing_ncurrindex = traindataNotMissing[ncurrindex]; 
	    
	          //if the mapping string is less than the number of data sets then 
	          //there are leading 0's will set for leading 0's not missing and absent
	          int numch = szmapping.length();
	          int numleading0 = numdatasets - numch;
	          for (int nj = 0; nj < numleading0; nj++)
	          {
		     traindataObservedValues_ncurrindex[nj] = false;
		     traindataNotMissing_ncurrindex[nj] = true;
		  }

	          int nmappedindex = numleading0; //starting from the leading 0 position
	          for (int nj = 0; nj < numch; nj++)
	          {
	             char ch = szmapping.charAt(nj);

	             if (ch == '0')
	             {
		        traindataObservedValues_ncurrindex[nmappedindex] = false;
		        traindataNotMissing_ncurrindex[nmappedindex] = true;
		     }
	             else if (ch=='1')
	             {
		        traindataObservedValues_ncurrindex[nmappedindex] = true;
		        traindataNotMissing_ncurrindex[nmappedindex] = true;
		     }
	             else
	             {
		        //missing data
		        traindataObservedValues_ncurrindex[nmappedindex] = false;
		        traindataNotMissing_ncurrindex[nmappedindex] = false;
		     }
	             nmappedindex++;
		  }
	       }      

	       traindataObservedIndexPair[nseq] = new int[numtime_nseq_m1];
	       int[] traindataObservedIndexPair_nseq = traindataObservedIndexPair[nseq];
	       boolean[] currvals = traindataObservedValues[traindataObservedIndex[0]];

	       for (int nindex = 0; nindex < numtime_nseq_m1; nindex++)//traindataObservedIndex_nseq_m1
	       {
	          boolean[] nextvals = traindataObservedValues[traindataObservedIndex[nindex+1]];
	          StringBuffer sb = new StringBuffer();
	          for (int nmark = 0; nmark < numdatasets;nmark++)
	          {
	             if (currvals[nmark]&&nextvals[nmark])
		     {
		        sb.append("1");
		     }
	             else
	             {
		        sb.append("0");
		     }
		  }

	          BigInteger theBigInteger = new BigInteger(sb.toString(),2);
	          Object obj = hmObserved.get(theBigInteger);
	          int ncurrobserved;
	          if (obj == null)
	          {
		     //first time we saw the combination storing index
		     hmObserved.put(theBigInteger,Integer.valueOf(nobserved));
	             ncurrobserved = nobserved;
	             nobserved++;

	             boolean[] pairflags = new boolean[numdatasets];
	             for (int nmark = 0; nmark < numdatasets;nmark++)
	             {
		        pairflags[nmark] = (currvals[nmark]&&nextvals[nmark]);
		     }
	             alobservedpairflags.add(pairflags);
		  }
	          else
	          {
		     //we already have seen this
		     ncurrobserved= ((Integer) hmObserved.get(theBigInteger)).intValue();
		  }
	          //storing which element pair the observation corresponds to
	          traindataObservedIndexPair_nseq[nindex] = ncurrobserved;
	          currvals = nextvals;
	       }
	    }//include
	}

	int numels = alobservedpairflags.size();

	//computes a tally for each flag combination observed of how frequently observed
	//int[] tallys = new int[numels];
	long[] tallys = new long[numels];

	//stores the total number of flag combinations observed
        //int ntotaltally = 0;
        long ntotaltally = 0;

	for (int nseq = 0; nseq < traindataObservedIndexPair.length; nseq++)
	{
	    if (bincludeseq[nseq])
	    {
               int[] traindataObservedIndexPair_nseq = traindataObservedIndexPair[nseq];
	       //updating element counts
	       for (int nindex = 0; nindex < traindataObservedIndexPair_nseq.length; nindex++)
	       {
	          tallys[traindataObservedIndexPair_nseq[nindex]]++;
	       }
	       //updates the total element counts
	       ntotaltally += traindataObservedIndexPair_nseq.length;
	    }
	}
       
	//first state always smoothed zero
	for (int nj = 0; nj < numdatasets; nj++)
        {
	   emissionprobs[0][nj][1] = dinformationsmooth*1.0/numbuckets;
           emissionprobs[0][nj][0] = 1- emissionprobs[0][nj][1];
	}
	

	//stores how many elements are currently assigned to the split node
	//int[] partitionTally = new int[numstates];
	long[] partitionTally = new long[numstates];

	//initally everything gets assigned to split node 0
	partitionTally[0] = ntotaltally;

	//stores the partition assignment of each split node
	// initially everything is assigned to split node 0
	int[] initStateAssign = new int[numels];

	//stores back pointers to the parent split nodes
	int[] backptr = new int[numstates];

	//stores for all prior splits how many elements would be split off
	//if partitioning on that mark
	//int[][] nextpartitionTally = new int[numstates-1][numdatasets];
	long[][] nextpartitionTally = new long[numstates-1][numdatasets];

	for (int niteration = 1; niteration < numstates; niteration++)
	{	
	    //need to the same number of assignments as we have states
 
	   for (int ni = 0; ni < niteration-1; ni++)
	   {
	       //initializes the nextpartitionTally to 0

	       //int[] nextpartitionTally_ni = nextpartitionTally[ni];
	       long[] nextpartitionTally_ni = nextpartitionTally[ni];
	       for (int nmark = 0; nmark < numdatasets; nmark++)
	       {
	          nextpartitionTally_ni[nmark] = 0;
	       }
	   }

	   for (int nel = 0; nel < tallys.length; nel++)
	   {
	       //considering each mark to partition on

	      boolean[] pairflags = (boolean[]) alobservedpairflags.get(nel);
	      //gets the current assignment of the element type
	      int initStateAssign_nel = initStateAssign[nel];

	      //counting how many from current partition would be split
	      for (int nsplitmark = 0; nsplitmark <numdatasets; nsplitmark++)
	      {
	         if (pairflags[nsplitmark])
		 {
		     //if element is positive for the mark then increments the count
		    nextpartitionTally[initStateAssign_nel][nsplitmark] += tallys[nel];
		 }
	      }
	   }

	   //finding the split that results in the greatest information change
	   double dbestinformationchange = 0;
	   //stores the best split
	   int nbestsplit = -1;
	   //stores the best split mark
	   int nbestsplitmark = -1;

           for (int nsplit = 0; nsplit < niteration; nsplit++)
           {
	       //considering all previous splits to split again

	       //fraction of total weighted elements in this node about to be split
	       //int partitionTally_nsplit = partitionTally[nsplit];
	       long partitionTally_nsplit = partitionTally[nsplit];
	       double dprobfull = partitionTally_nsplit/(double) ntotaltally;
	       double dprobfullterm;
	       if (dprobfull> 0)
	       {
		   //the information term for this node about to be split
                  dprobfullterm = dprobfull*Math.log(dprobfull);
	       }
	       else
	       {
		   dprobfullterm = 0;
	       }

	       //int[] nextpartitionTally_nsplit = nextpartitionTally[nsplit];
	       long[] nextpartitionTally_nsplit = nextpartitionTally[nsplit];

	       for (int nsplitmark = 0; nsplitmark < numdatasets; nsplitmark++)
	       {
		   //considering each mark to split on

		   //numerator is how many weighted elements would remain in the partition after splitting on the mark
	          double dprob1 = (partitionTally_nsplit-nextpartitionTally_nsplit[nsplitmark])/(double) ntotaltally;

		  //numerator is how many weighted elements would go to the new partition
	          double dprob2 = nextpartitionTally_nsplit[nsplitmark]/(double) ntotaltally;
	          double dinformationchange = dprobfullterm;
		  if (dprob1 > 0)
		  {
		      dinformationchange -= dprob1*Math.log(dprob1); 
		  }

		  if (dprob2 > 0)
		  {
		      dinformationchange -= dprob2*Math.log(dprob2); 
		  }
		  //-p_1*log(p_1)-p_2*log(p_2)- -(p1+p2)*log(p1+p2)

		  if (dinformationchange > dbestinformationchange)
		  {
	             dbestinformationchange = dinformationchange;
	             nbestsplit = nsplit;
	             nbestsplitmark = nsplitmark;
	          }
	       }
	   }

	    if (ChromHMM.BVERBOSE)
	    {
	       System.out.println("====>\t"+nbestsplit+"\t"+nbestsplitmark+"\t"+datasets[nbestsplitmark]+"\t"+dbestinformationchange+"\t"+
                                 nextpartitionTally[nbestsplit][nbestsplitmark]+"\t"+partitionTally[nbestsplit]);
	    }

	    if (nbestsplit == -1)
	    {
		throw new IllegalArgumentException("On this data the INFORMATION initialization strategy can only support "+niteration+" states. "+
                                                   "Check if the binarization was done correctly, and if so "+
						   "use the RANDOM or LOAD options for more states");
	    }

	    //the number of elements in the new split
	    //int numnewsplit =nextpartitionTally[nbestsplit][nbestsplitmark];
	    long numnewsplit =nextpartitionTally[nbestsplit][nbestsplitmark];
	    partitionTally[niteration] = numnewsplit;

	    //removes from the node we splitting from those we just split;
	    partitionTally[nbestsplit] -= numnewsplit;

	    //stores a back pointer to the node split to generate this one
	    backptr[niteration] = nbestsplit;

	    for (int nel = 0; nel < tallys.length; nel++)
	    {
		//goes through all elements and if recorded as being part of this split and positive for the split mark
		//then updates its initial state
	       boolean[] pairflags = (boolean[]) alobservedpairflags.get(nel);
	       if ((initStateAssign[nel]==nbestsplit)&& (pairflags[nbestsplitmark]))
	       {
		   initStateAssign[nel]= niteration;
	       }
	    }
	}

	//int[][] postally = new int[numstates][numdatasets];
       long[][] postally = new long[numstates][numdatasets];

       for (int nel = 0; nel < tallys.length; nel++)
       {
          boolean[] pairflags = (boolean[]) alobservedpairflags.get(nel);

	  for (int nmark = 0; nmark < numdatasets; nmark++)
	  {
	     if(pairflags[nmark])		
	     {
		 //this have a positive assignment for this flag
		 //walking to all ancestors and incrementing tally for this mark
	        int ncurrstate = initStateAssign[nel];
                do
	        {
		   postally[ncurrstate][nmark] += tallys[nel];		  
		   ncurrstate = backptr[ncurrstate];
		}
	        while (ncurrstate != 0);
	     }
	  }
       }

       //int[] partitionTallySum = new int[partitionTally.length];
       long[] partitionTallySum = new long[partitionTally.length];
       for (int nstate = 1; nstate < numstates; nstate++)
       {
	   //figures out how many descendants there are of each split
          int ncurrstate = nstate;
          //int ntotsum = 0;
          //int ncurrval = partitionTally[ncurrstate];
          long ncurrval = partitionTally[ncurrstate];
          do
          {
	      //incrementing counts for the parents
	     partitionTallySum[ncurrstate] += ncurrval;
      	     ncurrstate = backptr[ncurrstate];	    
	  }
       	  while (ncurrstate != 0);
       }

       for (int nstate = 1; nstate < numstates; nstate++)
       {
          for (int nj = 0; nj < numdatasets; nj++)
	  {
	     emissionprobs[nstate][nj][1] = dinformationsmooth*1.0/numbuckets
                                          +(1-dinformationsmooth)*(postally[nstate][nj]/(double) partitionTallySum[nstate]);
	     emissionprobs[nstate][nj][0] = 1- emissionprobs[nstate][nj][1];
	  }
       }

       //initialize the inital probability based on the partition of the first vector

       //int[] numstarts = new int[numstates];
       long[] numstarts = new long[numstates];
       //int nvalidseq = 0;
       long nvalidseq = 0;

       for (int nseq = 0; nseq < traindataObservedIndexPair.length; nseq++)
       {
	  if ((bincludeseq[nseq])&&(traindataObservedIndexPair[nseq].length > 0))
	  {
             //updated in v1.14
             numstarts[initStateAssign[traindataObservedIndexPair[nseq][0]]]++;
	     nvalidseq++;
	  }
       }
       
       for (int ni = 0; ni < probinit.length; ni++)
       {
	   //weighted probability of uniform and the fraction of starts from that partition
	   probinit[ni] = dinformationsmooth*1.0/numstates+(1-dinformationsmooth)*numstarts[ni]/(double) nvalidseq;  //updated in v1.14
       }


       //determining initial settings for the transition probability
       //int[][] transitiontally = new int[numstates][numstates];
       long[][] transitiontally = new long[numstates][numstates];
       int nnextstate;
       for (int nseq = 0; nseq < traindataObservedIndexPair.length; nseq++)
       {
	  if ((bincludeseq[nseq])&&(traindataObservedIndexPair[nseq].length > 0))
	  {
	      //going through each sequence
             int[] traindataObservedIndexPair_nseq = traindataObservedIndexPair[nseq];
	     int nprevstate = initStateAssign[traindataObservedIndexPair_nseq[0]];
             for (int nindex = 1; nindex < traindataObservedIndexPair_nseq.length; nindex++)
             {
	        //going through all the elements of the sequence computing a tallly of each bigram
	        nnextstate = initStateAssign[traindataObservedIndexPair_nseq[nindex]];
	        transitiontally[nprevstate][nnextstate]++;
	        nprevstate = nnextstate;
	     }
	  }
       }

       for (int ni = 0; ni < numstates; ni++)
       {
          double[] transitionprobs_ni = transitionprobs[ni];
	  double dnumfromi = 0;
	  //int[] transitiontally_ni = transitiontally[ni];
	  long[] transitiontally_ni = transitiontally[ni];

	  for (int nj = 0; nj < numstates; nj++)
	  {
	      dnumfromi += transitiontally_ni[nj];
	  }
          for (int nj = 0; nj < numstates; nj++)
          {
	     elim[ni][nj] = false;
	     //computes a weighted average of a uniform transition probability and the empirical frequency
	     transitionprobs_ni[nj] = dinformationsmooth*1.0/numstates
		                       +(1-dinformationsmooth)*(transitiontally_ni[nj]/(double) dnumfromi);
	     transitionprobsindex[ni][nj] = nj;
	     transitionprobsindexCol[ni][nj] = nj;
	  }
	  transitionprobsnum[ni] = numstates;
          transitionprobsnumCol[ni] = numstates;
       }
    }


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    public void informationInitializeNested() 
    {
       //inital probability vector
       probinit = new double[numstates];

       //creates the emission probability matrix
       emissionprobs = new double[numstates][numdatasets][numbuckets];

       //set to true if a transition has been eliminated
       elim = new boolean[numstates][numstates];

       //initalize the transition matrix
       transitionprobs = new double[numstates][numstates];

       //initialize index of the next non-zero transition
       transitionprobsindex = new int[numstates][numstates];

       //initalize number of non-zero transitions
       transitionprobsnum = new int[numstates];

       //initalize column-wise index of non-zero transitions
       transitionprobsindexCol = new int[numstates][numstates];
 
       //number of non-zero column transitions
       transitionprobsnumCol = new int[numstates];

       int[][] traindataObservedIndexPair = new int[traindataObservedIndex.length][];
	
       ArrayList alobservedpairflags = new ArrayList(); //is an index from the element combination to the associated flags

       HashMap hmObserved= new HashMap(); //is an index from the associated flags to the element index

       int nobserved= 0; //index on the unique observation combination we are observing

       int[] samples= null;

       int ncurrnumincludeseq;
       if (numincludeseq >= 1)
       {
	   ncurrnumincludeseq = Math.min(numincludeseq, chromfiles.length);
	   samples = new int[ncurrnumincludeseq];
       }
       else
       {
	   ncurrnumincludeseq = chromfiles.length;
       }

       boolean[] bincludeseq = new boolean[chromfiles.length];
       for (int nk = 0; nk < bincludeseq.length; nk++)
       {
	   bincludeseq[nk] = true;
       }

       if (numincludeseq >=1)
       {
          for (int ni = 0; ni < samples.length; ni++)
	  {
	     samples[ni] = ni;
	  }
	  //Random theRandom2 = new Random(412);//433

	  //sampling without replacement data locations wanted
	  for (int ni = ncurrnumincludeseq; ni < chromfiles.length; ni++)
          {
	     if (theRandom.nextDouble() < ncurrnumincludeseq/((double) ni +1))
	     {
      	        samples[theRandom.nextInt(ncurrnumincludeseq)] = ni;
      	     }
	  }

          for (int nk = 0; nk < bincludeseq.length; nk++)
          {
             bincludeseq[nk] = false;
	  }

          for (int nk = 0; nk < samples.length; nk++)
	  {
             bincludeseq[samples[nk]] = true;
		   //System.out.println("including\t"+samples[nk]+"\t"+bincludeseq[samples[nk]]);
          }
       }

       //generates all vectors of combinations of consecutive of '1' calls
       for (int nseq = 0; nseq <traindataObservedIndex.length; nseq++)
       {
          if (bincludeseq[nseq])
	  {
	     //going through each sequence

	     int[] traindataObservedIndex_nseq = traindataObservedIndex[nseq];
	     int traindataObservedIndex_nseq_m1 = traindataObservedIndex_nseq.length -1;
	     traindataObservedIndexPair[nseq] = new int[traindataObservedIndex_nseq_m1];
	     int[] traindataObservedIndexPair_nseq = traindataObservedIndexPair[nseq];
	     boolean[] currvals = traindataObservedValues[traindataObservedIndex_nseq[0]];

	     for (int nindex = 0; nindex <  traindataObservedIndex_nseq_m1; nindex++)
	     {
	        boolean[] nextvals = traindataObservedValues[traindataObservedIndex_nseq[nindex+1]];
	        StringBuffer sb = new StringBuffer();
	        for (int nmark = 0; nmark < numdatasets;nmark++)
	        {
	           if (currvals[nmark]&&nextvals[nmark])
		   {
		      sb.append("1");
		   }
	           else
	           {
		      sb.append("0");
		   }
	        }

	        BigInteger theBigInteger = new BigInteger(sb.toString(),2);
	        Object obj = hmObserved.get(theBigInteger);
	        int ncurrobserved;
	        if (obj == null)
	        {
		  //first time we saw the combination storing index
		   hmObserved.put(theBigInteger,Integer.valueOf(nobserved));
	           ncurrobserved = nobserved;
	           nobserved++;

	           boolean[] pairflags = new boolean[numdatasets];
	           for (int nmark = 0; nmark < numdatasets;nmark++)
	           {
		      pairflags[nmark] = (currvals[nmark]&&nextvals[nmark]);
		   }
	           alobservedpairflags.add(pairflags);
	        }
	        else
	        {
		   //we already have seen this
		   ncurrobserved= ((Integer) hmObserved.get(theBigInteger)).intValue();
	        }

	        //storing which element pair the observation corresponds to
	        traindataObservedIndexPair_nseq[nindex] = ncurrobserved;
	        currvals = nextvals;
	     }
	  }
       }

       int numels = alobservedpairflags.size();

       //computes a tally for each flag combination observed of how frequently observed
       //int[] tallys = new int[numels];
       long[] tallys = new long[numels];

       //stores the total number of flag combinations observed
       //int ntotaltally = 0;
       long ntotaltally = 0;

       for (int nseq = 0; nseq < traindataObservedIndexPair.length; nseq++)
       {
	    if (bincludeseq[nseq])
	    {
               int[] traindataObservedIndexPair_nseq = traindataObservedIndexPair[nseq];
	       //updating element counts
	       for (int nindex = 0; nindex < traindataObservedIndexPair_nseq.length; nindex++)
	       {
	          tallys[traindataObservedIndexPair_nseq[nindex]]++;
	       }
	       //updates the total element counts
	       ntotaltally += traindataObservedIndexPair_nseq.length;
	    }
	}
       
	//first state always smoothed zero
	for (int nj = 0; nj < numdatasets; nj++)
        {
	   emissionprobs[0][nj][1] = dinformationsmooth*1.0/numbuckets;
           emissionprobs[0][nj][0] = 1- emissionprobs[0][nj][1];
	}
	

	//stores how many elements are currently assigned to the split node
	//int[] partitionTally = new int[numstates];
	long[] partitionTally = new long[numstates];

	//initally everything gets assigned to split node 0
	partitionTally[0] = ntotaltally;

	//stores the partition assignment of each split node
	// initially everything is assigned to split node 0
	int[] initStateAssign = new int[numels];

	//stores back pointers to the parent split nodes
	int[] backptr = new int[numstates];

	//stores for all prior splits how many elements would be split off
	//if partitioning on that mark
	//int[][] nextpartitionTally = new int[numstates-1][numdatasets];
	long[][] nextpartitionTally = new long[numstates-1][numdatasets];

	for (int niteration = 1; niteration < numstates; niteration++)
	{	
	    //need to the same number of assignments as we have states
 
	   for (int ni = 0; ni < niteration-1; ni++)
	   {
	       //initializes the nextpartitionTally to 0

	       //int[] nextpartitionTally_ni = nextpartitionTally[ni];
	       long[] nextpartitionTally_ni = nextpartitionTally[ni];
	       for (int nmark = 0; nmark < numdatasets; nmark++)
	       {
	          nextpartitionTally_ni[nmark] = 0;
	       }
	   }

	   for (int nel = 0; nel < tallys.length; nel++)
	   {
	       //considering each mark to partition on

	      boolean[] pairflags = (boolean[]) alobservedpairflags.get(nel);
	      //gets the current assignment of the element type
	      int initStateAssign_nel = initStateAssign[nel];

	      //counting how many from current partition would be split
	      for (int nsplitmark = 0; nsplitmark <numdatasets; nsplitmark++)
	      {
	         if (pairflags[nsplitmark])
		 {
		     //if element is positive for the mark then increments the count
		    nextpartitionTally[initStateAssign_nel][nsplitmark] += tallys[nel];
		 }
	      }
	   }

	   //finding the split that results in the greatest information change
	   double dbestinformationchange = 0;
	   //stores the best split
	   int nbestsplit = -1;
	   //stores the best split mark
	   int nbestsplitmark = -1;

           for (int nsplit = 0; nsplit < niteration; nsplit++)
           {
	       //considering all previous splits to split again

	       //fraction of total weighted elements in this node about to be split
	       //int partitionTally_nsplit = partitionTally[nsplit];
	       long partitionTally_nsplit = partitionTally[nsplit];
	       double dprobfull = partitionTally_nsplit/(double) ntotaltally;
	       double dprobfullterm;
	       if (dprobfull> 0)
	       {
		   //the information term for this node about to be split
                  dprobfullterm = dprobfull*Math.log(dprobfull);
	       }
	       else
	       {
		   dprobfullterm = 0;
	       }

	       //int[] nextpartitionTally_nsplit = nextpartitionTally[nsplit];
	       long[] nextpartitionTally_nsplit = nextpartitionTally[nsplit];

	       for (int nsplitmark = 0; nsplitmark < numdatasets; nsplitmark++)
	       {
		   //considering each mark to split on

		   //numerator is how many weighted elements would remain in the partition after splitting on the mark
	          double dprob1 = (partitionTally_nsplit-nextpartitionTally_nsplit[nsplitmark])/(double) ntotaltally;

		  //numerator is how many weighted elements would go to the new partition
	          double dprob2 = nextpartitionTally_nsplit[nsplitmark]/(double) ntotaltally;
	          double dinformationchange = dprobfullterm;
		  if (dprob1 > 0)
		  {
		      dinformationchange -= dprob1*Math.log(dprob1); 
		  }

		  if (dprob2 > 0)
		  {
		      dinformationchange -= dprob2*Math.log(dprob2); 
		  }
		  //-p_1*log(p_1)-p_2*log(p_2)- -(p1+p2)*log(p1+p2)

		  if (dinformationchange > dbestinformationchange)
		  {
	             dbestinformationchange = dinformationchange;
	             nbestsplit = nsplit;
	             nbestsplitmark = nsplitmark;
	          }
	       }
	    }

	    if (ChromHMM.BVERBOSE)
	    {
	       System.out.println("====>\t"+nbestsplit+"\t"+nbestsplitmark+"\t"+datasets[nbestsplitmark]+"\t"+dbestinformationchange+"\t"+
                                 nextpartitionTally[nbestsplit][nbestsplitmark]+"\t"+partitionTally[nbestsplit]);
	    }

	    if (nbestsplit == -1)
	    {
		throw new IllegalArgumentException("On this data the INFORMATION initialization strategy can only support "+niteration+" states. "+
                                                   "Check if the binarization was done correctly, and if so "+
						   "use the RANDOM or LOAD options for more states");
	    }

	    //the number of elements in the new split
	    //int numnewsplit =nextpartitionTally[nbestsplit][nbestsplitmark];
	    long numnewsplit = nextpartitionTally[nbestsplit][nbestsplitmark];
	    partitionTally[niteration] = numnewsplit;

	    //removes from the node we splitting from those we just split;
	    partitionTally[nbestsplit] -= numnewsplit;

	    //stores a back pointer to the node split to generate this one
	    backptr[niteration] = nbestsplit;

	    for (int nel = 0; nel < tallys.length; nel++)
	    {
		//goes through all elements and if recorded as being part of this split and positive for the split mark
		//then updates its initial state
	       boolean[] pairflags = (boolean[]) alobservedpairflags.get(nel);
	       if ((initStateAssign[nel]==nbestsplit)&& (pairflags[nbestsplitmark]))
	       {
		   initStateAssign[nel]= niteration;
	       }
	    }
	}

       //int[][] postally = new int[numstates][numdatasets];
       long[][] postally = new long[numstates][numdatasets];

       for (int nel = 0; nel < tallys.length; nel++)
       {
          boolean[] pairflags = (boolean[]) alobservedpairflags.get(nel);

	  for (int nmark = 0; nmark < numdatasets; nmark++)
	  {
	     if(pairflags[nmark])		
	     {
		 //this have a positive assignment for this flag
		 //walking to all ancestors and incrementing tally for this mark
	        int ncurrstate = initStateAssign[nel];
                do
	        {
		   postally[ncurrstate][nmark] += tallys[nel];		  
		   ncurrstate = backptr[ncurrstate];
		}
	        while (ncurrstate != 0);
	     }
	  }
       }

       //int[] partitionTallySum = new int[partitionTally.length];
       long[] partitionTallySum = new long[partitionTally.length];
       for (int nstate = 1; nstate < numstates; nstate++)
       {
	   //figures out how many descendants there are of each split
          int ncurrstate = nstate;
          //int ntotsum = 0;
          //int ncurrval = partitionTally[ncurrstate];
          long ncurrval = partitionTally[ncurrstate];
          do
          {
	      //incrementing counts for the parents
	     partitionTallySum[ncurrstate] += ncurrval;
      	     ncurrstate = backptr[ncurrstate];	    
	  }
       	  while (ncurrstate != 0);
       }

       for (int nstate = 1; nstate < numstates; nstate++)
       {
          for (int nj = 0; nj < numdatasets; nj++)
	  {
	     emissionprobs[nstate][nj][1] = dinformationsmooth*1.0/numbuckets
                                          +(1-dinformationsmooth)*(postally[nstate][nj]/(double) partitionTallySum[nstate]);
	     emissionprobs[nstate][nj][0] = 1- emissionprobs[nstate][nj][1];
	  }
       }

       //initialize the inital probability based on the partition of the first vector

       //int[] numstarts = new int[numstates];
       long[] numstarts = new long[numstates];
       //int nvalidseq = 0;
       long nvalidseq = 0;
       for (int nseq = 0; nseq < traindataObservedIndexPair.length; nseq++)
       {
	  if ((bincludeseq[nseq])&&(traindataObservedIndexPair[nseq].length > 0))
	  {
	      //updated in v1.14
             numstarts[initStateAssign[traindataObservedIndexPair[nseq][0]]]++;
	     nvalidseq++; 
	  }
       }
       
       for (int ni = 0; ni < probinit.length; ni++)
       {
	   //weighted probability of uniform and the fraction of starts from that partition
	   probinit[ni] = dinformationsmooth*1.0/numstates+(1-dinformationsmooth)*numstarts[ni]/(double) nvalidseq; //updated in v1.14
       }


       //determining initial settings for the transition probability
       //int[][] transitiontally = new int[numstates][numstates];
       long[][] transitiontally = new long[numstates][numstates];
       int nnextstate;
       for (int nseq = 0; nseq < traindataObservedIndexPair.length; nseq++)
       {
	  if ((bincludeseq[nseq])&&(traindataObservedIndexPair[nseq].length > 0)) //updated in v1.14
	  {
	     //going through each sequence
             int[] traindataObservedIndexPair_nseq = traindataObservedIndexPair[nseq];
	     int nprevstate = initStateAssign[traindataObservedIndexPair_nseq[0]];
             for (int nindex = 1; nindex < traindataObservedIndexPair_nseq.length; nindex++)
             {
	         //going through all the elements of the sequence computing a tallly of each bigram
	        nnextstate = initStateAssign[traindataObservedIndexPair_nseq[nindex]];
	        transitiontally[nprevstate][nnextstate]++;
	        nprevstate = nnextstate;
	     }
	  }
       }

       for (int ni = 0; ni < numstates; ni++)
       {
          double[] transitionprobs_ni = transitionprobs[ni];
	  double dnumfromi = 0;
	  //int[] transitiontally_ni = transitiontally[ni];
	  long[] transitiontally_ni = transitiontally[ni];

	  for (int nj = 0; nj < numstates; nj++)
	  {
	      dnumfromi += transitiontally_ni[nj];
	  }
          for (int nj = 0; nj < numstates; nj++)
          {
	     elim[ni][nj] = false;
	     //computes a weighted average of a uniform transition probability and the empirical frequency
	     transitionprobs_ni[nj] = dinformationsmooth*1.0/numstates
		                       +(1-dinformationsmooth)*(transitiontally_ni[nj]/(double) dnumfromi);
	     transitionprobsindex[ni][nj] = nj;
	     transitionprobsindexCol[ni][nj] = nj;
	  }
	  transitionprobsnum[ni] = numstates;
          transitionprobsnumCol[ni] = numstates;
       }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////


    /**
     * Randomly initializes parameters of the HMM based on a uniform distribution
     */
    public void randomlyInitializeParams()
    {
	//initializes the initial state parameters
	double dsum = 0;
	probinit = new double[numstates];
	for (int ni = 0; ni < probinit.length; ni++)
	{
	    probinit[ni] = theRandom.nextDouble();
	    dsum += probinit[ni];
	}

	for (int ni = 0; ni < probinit.length; ni++)
	{
	    probinit[ni] /= dsum;
	}

	//set to true if a transition has been eliminated
	elim = new boolean[numstates][numstates];

	//initalize the transition matrix
	transitionprobs = new double[numstates][numstates];

	//initialize index of the next non-zero transition
	transitionprobsindex = new int[numstates][numstates];

	//initalize number of non-zero transitions
	transitionprobsnum = new int[numstates];

	//initalize column-wise index of non-zero transitions
        transitionprobsindexCol = new int[numstates][numstates];

	//number of non-zero column transitions
        transitionprobsnumCol = new int[numstates];

	//uniformly randomly assigns transition probability values
	//also initializes transition indicies
	for (int ni = 0; ni < transitionprobs.length; ni++)
	{
	    dsum = 0;
	    double[] transitionprobs_ni = transitionprobs[ni];
	    for (int nj = 0; nj < transitionprobs_ni.length; nj++)
	    {
		elim[ni][nj] = false;
		double dval = theRandom.nextDouble();
		transitionprobs_ni[nj] = dval;
		dsum += transitionprobs_ni[nj];
		transitionprobsindex[ni][nj] = nj;
		transitionprobsindexCol[ni][nj] = nj;
	    }
	    transitionprobsnum[ni] = numstates;
	    transitionprobsnumCol[ni] = numstates;

	    for (int nj = 0; nj < transitionprobs_ni.length; nj++)
	    {
		transitionprobs_ni[nj] /= dsum;
	    }
	}

	//uniformly randomly assigns emission probability values
	//also initializes emission indicies
	emissionprobs = new double[numstates][numdatasets][numbuckets];
	for (int ni = 0; ni < emissionprobs.length; ni++)
	{
	    double[][] emissionprobs_ni = emissionprobs[ni];
	    for (int nj = 0; nj < emissionprobs_ni.length; nj++)
	    {
		double[] emissionprobs_ni_nj = emissionprobs_ni[nj];
		dsum = 0;
		for (int nk = 0; nk < emissionprobs_ni_nj.length; nk++)
		{
		    double dval = theRandom.nextDouble();
		    dsum += dval;
		    emissionprobs_ni_nj[nk] = dval;
		}

		for (int nk = 0; nk < emissionprobs_ni_nj.length; nk++)
		{
		    emissionprobs_ni_nj[nk] /= dsum;
		}
	    }
	}
    }


    ///////////////////////////////////////////////////////////////////////////////////////////

    /**
     * loadModel loads a model directly by calling loadModelSmooth with emission parameters set to 0
     */
    public void loadModel() throws IOException
    {
	loadModelSmooth(0,0);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////

    /**
     * Loads the model contained in szInitFile with the smoothing specified
     * based on the dsmoothtransition and dsmoothemission parameters
     */
    public void loadModelSmooth(double dproceduresmoothemission, double dproceduresmoothtransition) throws IOException
    {
	BufferedReader br =  Util.getBufferedReader(szInitFile);
	String szLine;
	szLoadHeader = br.readLine();
	if (szLoadHeader == null)
	{
	    throw new IllegalArgumentException(szLoadHeader+" is empty!");
	}
        StringTokenizer st = new StringTokenizer(szLoadHeader,"\t");

	//first token of the first line of the model file is assume to give the number of states
	numstates = Integer.parseInt(st.nextToken().trim());
	numdatasets = Integer.parseInt(st.nextToken().trim());
	if (datasets ==null)
	{
	    datasets = new String[numdatasets];
	}
	chorder = st.nextToken().trim().charAt(0);

	if ((nstateorder != ChromHMM.STATEORDER_TRANSITION)&&(nstateorder != ChromHMM.STATEORDER_EMISSION)&&(nstateorder != ChromHMM.STATEORDER_FIXED))
	{
	   nstateorder = -1;
	   for (int ni = 0; ni < ChromHMM.ORDERCHARS.length; ni++)
	   {
	      if (chorder == ChromHMM.ORDERCHARS[ni])
	      {
	         nstateorder = ni;
	         break;
	      }
	   }

	   if (nstateorder == -1)
	   {
	      throw new IllegalArgumentException(chorder+" is an invalid order type");
	   }
	   szorder = ChromHMM.ORDERSTRINGS[nstateorder];
	}

	probinit = new double[numstates];

	//transitions that are strictly 0 are elminated
        elim = new boolean[numstates][numstates];

	transitionprobs = new double[numstates][numstates];

	//contains the indicies of non-0 transitions
        transitionprobsindex = new int[numstates][numstates];
        transitionprobsnum = new int[numstates];

        transitionprobsindexCol = new int[numstates][numstates];
        transitionprobsnumCol = new int[numstates];

	//initialize the transition structure without any parameter elimination
        for (int ni = 0; ni < transitionprobs.length; ni++)
        {
	   for (int nj = 0; nj < transitionprobs[ni].length; nj++)
           {
	      elim[ni][nj] = false;
	      transitionprobsindex[ni][nj] = nj;
	      transitionprobsindexCol[ni][nj] = nj;
           }
           transitionprobsnum[ni] = numstates;
	   transitionprobsnumCol[ni] = numstates;
	}

	emissionprobs = new double[numstates][numdatasets][numbuckets];

	boolean btransition0 = false;

	while ((szLine = br.readLine())!=null)
	{
	    st = new StringTokenizer(szLine,"\t");
	    String szvartype = st.nextToken().trim();

	    if (szvartype.equalsIgnoreCase("probinit"))
	    {
		//reading an inital probability
		int nstate = Integer.parseInt(st.nextToken().trim())-1;
		double dprob = Double.parseDouble(st.nextToken().trim());
		probinit[nstate] = dprob;
	    }
	    else if (szvartype.equalsIgnoreCase("transitionprobs"))
	    {
		int nfrom = Integer.parseInt(st.nextToken().trim())-1;
		int nto = Integer.parseInt(st.nextToken().trim())-1;
		double dprob = Double.parseDouble(st.nextToken().trim());
		//this smooths the transition probability if dproceduresmmothtransition>0 using a weighted average with uniform 
		transitionprobs[nfrom][nto] = dproceduresmoothtransition/((double) transitionprobs.length)+(1-dproceduresmoothtransition)*dprob;
		if (transitionprobs[nfrom][nto] == 0)
		{
		    //we have a 0 transition
		    btransition0 = true;
		    elim[nfrom][nto] = true;
		}
	    }
	    else if (szvartype.equalsIgnoreCase("emissionprobs"))
	    {
		int nstate = Integer.parseInt(st.nextToken().trim())-1;
		int nmod = Integer.parseInt(st.nextToken().trim());
		String szdataset = st.nextToken().trim();
		if (datasets[nmod]==null)
		{
		    datasets[nmod] = szdataset;//st.nextToken().trim();
		}
		else if (!szdataset.equals(datasets[nmod]))
		{
		    throw new IllegalArgumentException("For mark "+nmod+" in model file found "+szdataset+", but expecting "+datasets[nmod]+"."
                                                      +" This might be because columns of binarized file are in a different order from model learning.");
		}
		int nval = Integer.parseInt(st.nextToken().trim());
		double dprob = Double.parseDouble(st.nextToken().trim());
		//smooths the emission probability if dproceduresmoothemission>0 using a weighted average with uniform
		emissionprobs[nstate][nmod][nval] = dproceduresmoothemission/((double) numbuckets)+(1-dproceduresmoothemission)*dprob;
	    }
	    else
	    {
		throw new IllegalArgumentException(szvartype+" is not recognized in the input model file");
	    }
	}
	br.close();

	if (btransition0)
        {
	   //we have a non-zero transition will update the sparse indicies
	   for (int ni = 0; ni < transitionprobs.length; ni++)
	   {
	       int nindex = 0;		       
	       boolean[] elim_ni = elim[ni];
	       int[] transitionprobsindex_ni = transitionprobsindex[ni];
	       for (int nj = 0; nj < transitionprobsindex_ni.length; nj++)
	       {
	          if (!elim_ni[nj])
		  {
		     //we have not eliminated this transition
	             //stores its index in order and add sum to denominator
		     transitionprobsindex_ni[nindex] = nj;
	       	     nindex++;
		  }
	       }
	       //update the number of valid transitions
	       transitionprobsnum[ni] = nindex; 
	    }

	    for (int ni = 0; ni < transitionprobs.length; ni++)
	    {
	       int nindex =0;
	       for (int nj = 0; nj < numstates; nj++)
	       {
	          if (!elim[nj][ni])
	          {
		     //copy into the column i the indicies of all non-eliminated transitions into i
		     transitionprobsindexCol[ni][nindex] = nj;
	             nindex++;
		  }
	       }
	       //updates the number of non-zero transitions from column i
	       transitionprobsnumCol[ni] = nindex; 
	    }
	}
    }

    /**
     * Takes an existing model and segmentation and outputs a confusion matrix by a selected subset
     */
    public void makeSegmentationConfusion() throws IOException
    {
        NumberFormat nf = NumberFormat.getInstance();
        nf.setMaximumFractionDigits(4);

       //number of non-zero transition required to be less than this at the more stringent cutoff 
       //for trying to exploit sparsity in the transition matrix for efficiency gains
       int nsparsecutoff = (int) (numstates * ChromHMM.SPARSECUTOFFRATIO);

       int[] numtime = new int[traindataObservedIndex.length];

       //stores the maximum number of locations in any sequence and in each sequence
       int nmaxtime = 0;
       for (int nseq = 0; nseq < traindataObservedIndex.length; nseq++)
       {
          numtime[nseq] = traindataObservedIndex[nseq].length;
          if (numtime[nseq] > nmaxtime)
	  {
      	     nmaxtime = numtime[nseq];
	  }
       }

       //double
       double[][] fullposterior = null;

       int[] fullmax = null;
       double[][] confusion = new double[numstates][numstates];
       double[][] normalizedconfusion = new double[numstates][numstates];

       if (breadposterior)
       {
          fullposterior = new double[nmaxtime][numstates];
       }

       if ((breadstatebyline)||(breadsegment))
       {
	  fullmax = new int[nmaxtime];
       }


       if (ChromHMM.BVERBOSE)
       {
          System.out.println("Maximum number of locations\t"+nmaxtime);
       }


       //stores the emission probability for the i^th combination of marks in the j^th state
       double[][] emissionproducts = new double[traindataObservedValues.length][numstates];

       //stores temporary product terms
       double[] tempproductbetaemiss = new double[numstates];

       //This stores the alpha values at each time point and number of states
       double[][] alpha = new double[nmaxtime][numstates];

       //Temporary storage of the gamma's for each state
       double[][] gamma = new double[nmaxtime][numstates];

       //Temporary storage of the beta values for each state
       double[] beta_nt = new double[numstates];

       //Temporary storage of the beta values for each state at the next time point
       double[] beta_ntp1 = new double[numstates];

       //stores the scaling value for each time point
       double[] scale = new double[nmaxtime];

       //stores the transition probabilities for each column
       double[][] coltransitionprobs = new double[numstates][numstates];

       boolean[] includemarks = new boolean[numdatasets];

       double[] surplus = new double[numstates];
       double[] deficit = new double[numstates];
       double[] dstatesagree = new double[numstates];

       if (szincludemarks.length()!=numdatasets)
       {
	   throw new IllegalArgumentException("Number of marks in "+szincludemarks+" of "+szincludemarks.length()+" does not equal expected number of "+numdatasets);
       }
      

       if ((breadstatebyline)||(breadsegment))
       {
	   //stores the maximum assignment with all marks
          fullmax = new int[nmaxtime];
       }
       else
       {
	   //stores the posterior assignment with all marks
          fullposterior = new double[nmaxtime][numstates];
       }

       String szdatasets = "";
       for (int nmark = 0; nmark < includemarks.length; nmark++)
       {
	   if (szincludemarks.charAt(nmark) == '1')
	   {
	       //stores in includemarks those data sets that have a '1' for the mark
	       includemarks[nmark] = true;
	       if (szdatasets.equals(""))
	       {
		   szdatasets += datasets[nmark];
	       }
	       else
	       {
		   szdatasets += "," + datasets[nmark];
	       }
	   }
	   else if (szincludemarks.charAt(nmark) == '0')
	   {
	       includemarks[nmark] = false;
	   }
	   else
	   {
	       throw new IllegalArgumentException(szincludemarks+" is not a valid bit string for includemarks!");
	   }
       }

       RecIntString[] ordered = new RecIntString[chromfiles.length];
       for (int nindex = 0; nindex < ordered.length; nindex++)
       {
	   ordered[nindex] = new RecIntString(nindex,chromfiles[nindex]);
       }
       Arrays.sort(ordered,new RecIntStringCompare());



       hsprefix = new HashSet();

       for (int nseq = 0; nseq < traindataObservedIndex.length; nseq++)
       {
          int nordered_nseq = ordered[nseq].nindex;
	  //goes through each sequence

          int[] traindataObservedIndex_nseq = traindataObservedIndex[nordered_nseq];
          boolean[] traindataObservedSeqFlags_nseq = traindataObservedSeqFlags[nordered_nseq];

	  String szprefix = "";
	  if (!cellSeq[nordered_nseq].equals(""))
	  {
	     szprefix += cellSeq[nordered_nseq]+"_";
	  }
	  szprefix += numstates;
	  if (!szoutfileID.equals(""))
	  {
	     szprefix += "_"+szoutfileID;
          }
	  hsprefix.add(szprefix);
	  
	  if (breadposterior)
	  {
	     BufferedReader brprobs = null;
	     //creates the posterior file
	     String szposteriorinfilename = szsegmentdir+"/POSTERIOR/"+szprefix+"_"+chromSeq[nordered_nseq]+ChromHMM.SZPOSTERIOREXTENSION;

	     //updated in v.1.11 to check for  gz extensions
	     File f = new File (szposteriorinfilename);
	     File fgz = new File (szposteriorinfilename+".gz");
	     if ((!f.exists())&&(fgz.exists()))
	     {
		 szposteriorinfilename = szposteriorinfilename + ".gz";
	     }
	     brprobs = Util.getBufferedReader(szposteriorinfilename);// new BufferedReader(new FileReader(szposteriorinfilename));

	     //skips the header lines
	     brprobs.readLine();
	     brprobs.readLine();
	     String szLinePosterior;

	     int nline = 0;
	     while ((szLinePosterior = brprobs.readLine())!=null)
	     {
		StringTokenizer stposterior = new StringTokenizer(szLinePosterior,"\t ");
	        for (int nstate = 0; nstate < numstates; nstate++)
	        {
		   fullposterior[nline][nstate] = Double.parseDouble(stposterior.nextToken());
		}
		nline++;
	     }
	     brprobs.close(); 
	  } 
          else if (breadstatebyline)
	  {
	      String szcurrchrom = chromSeq[nordered_nseq];
	      //reads a file which has the state with the maximum posterior probability
	      String szmaxinfilename = szsegmentdir+"/STATEBYLINE/"+szprefix+"_"+szcurrchrom+ChromHMM.SZSTATEBYLINEEXTENSION;

	      //updated in v.1.11 to check for  gz extensions
	      File f = new File (szmaxinfilename);
	      File fgz = new File (szmaxinfilename+".gz");
	      if ((!f.exists())&&(fgz.exists()))
	      {
	         szmaxinfilename = szmaxinfilename + ".gz";
	      }
 	      BufferedReader brmax = Util.getBufferedReader(szmaxinfilename);// new BufferedReader(new FileReader(szmaxinfilename));
	      //skip the header lines
	      brmax.readLine();
	      brmax.readLine();
	      String szLineMax;
	      int nline = 0;
	      while ((szLineMax = brmax.readLine())!=null)
	      {
		  fullmax[nline] = Integer.parseInt(szLineMax)-1;
		  nline++;
	      }
	      brmax.close(); 
	  }
	  else if (breadsegment)
	  {

	     BufferedReader brbed = null;
	     //creates a file which has the maximum segmentation
	     //we only have one file per cell type here
	     String szcurrchrom = chromSeq[nordered_nseq];

	     String szsegmentinfilename = szsegmentdir+"/" + szprefix+ChromHMM.SZSEGMENTEXTENSION;

	     //updated in v.1.11 to check for  gz extensions
	     File f = new File (szsegmentinfilename);
	     File fgz = new File (szsegmentinfilename+".gz");
	     if ((!f.exists())&&(fgz.exists()))
	     {
                szsegmentinfilename = szsegmentinfilename + ".gz";
             }

	     brbed = Util.getBufferedReader(szsegmentinfilename); // new BufferedReader(new FileReader(szsegmentinfilename));
		 
	     String szLineMax;
	     while ((szLineMax = brbed.readLine())!=null)
	     {
		 StringTokenizer stchrom = new StringTokenizer(szLineMax,"\t");
		 String szchrom = stchrom.nextToken().trim();

		 if (szchrom.equals(szcurrchrom))
	         {
		     int nbegin = Integer.parseInt(stchrom.nextToken().trim())/nbinsize;
		     int nend = (Integer.parseInt(stchrom.nextToken().trim())-1)/nbinsize;
		     int nstate = Integer.parseInt(stchrom.nextToken().trim().substring(1))-1;		        		
		     for (int nj = nbegin; nj <= nend; nj++)
		     {
			 fullmax[nj] = nstate;
		     }
		 }	      
	     }
	     brbed.close();		
	  }


          if (bscaleemissions)
	  {
	     for (int ni = 0; ni < emissionproducts.length; ni++)
	     {
	        //going through each combination of marks
		if (traindataObservedSeqFlags_nseq[ni])
		{
		   //this signature of marks is observed on the current chromosome so
		   //updating its emission probabilities
		   double[] emissionproducts_ni = emissionproducts[ni];
		   boolean[] traindataObservedValues_ni = traindataObservedValues[ni];
		   boolean[] traindataNotMissing_ni = traindataNotMissing[ni];

		   for (int ns = 0; ns < numstates; ns++)
		   {
		      emissionproducts_ni[ns] = 1;
		   }

		   for (int nmod = 0; nmod < numdatasets; nmod++)
		   {
		      for (int ns = 0; ns < numstates; ns++)
		      {
			 if ((traindataNotMissing_ni[nmod])&&(includemarks[nmod]))
			 {
			    //we are include this marks emission probability
			    if (traindataObservedValues_ni[nmod])
			    {
			       emissionproducts_ni[ns] *= emissionprobs[ns][nmod][1];
			    }
			    else
			    {
			       emissionproducts_ni[ns] *= emissionprobs[ns][nmod][0];
			    }
			 }
			 // otherwise treated as missing omitting from product
		      }

		      double dmaxval = 0;
		      for (int ns = 0; ns < numstates; ns++)
		      {
		         if (emissionproducts_ni[ns] > dmaxval)
			 {
			    dmaxval = emissionproducts_ni[ns];
			 }
		      }

		      if (dmaxval <= 0)
		      {
		         for (int ns = 0; ns < numstates; ns++)
		         {
		            emissionproducts_ni[ns] = 1;
		         }
		      }
		      else
		      {
		         for (int ns = 0; ns < numstates; ns++)
		         {
		            emissionproducts_ni[ns]/= dmaxval;
		         }
		      }
		   }
		}
	     }
	  }
	  else
	  {
	     for (int ni = 0; ni < emissionproducts.length; ni++)
             {
	        //going through each combination of marks
	        if (traindataObservedSeqFlags_nseq[ni])
	        {
	           //this signature of marks is observed on the current chromosome so
	           //updating its emission probabilities
	           double[] emissionproducts_ni = emissionproducts[ni];
	           boolean[] traindataObservedValues_ni = traindataObservedValues[ni];
	           boolean[] traindataNotMissing_ni = traindataNotMissing[ni];	
		   boolean ballzero = true;
	  
	           for (int ns = 0; ns < numstates; ns++)
	           {
	              double dproduct = 1;
	              double[][] emissionprobs_ni = emissionprobs[ns];

		      //going through all marks
		      for (int nmod = 0; nmod < numdatasets; nmod++)
	              {
		         if ((traindataNotMissing_ni[nmod])&&(includemarks[nmod]))
		         {
			     //we have observed the mark
		            if (traindataObservedValues_ni[nmod])
		            {
		               dproduct *= emissionprobs_ni[nmod][1];
		            }
		            else 
	                    {
		               dproduct *= emissionprobs_ni[nmod][0];
		            }
			 }
		         // otherwise treated as missing omitting from product
		      }
	              emissionproducts_ni[ns] = dproduct;

		      if (dproduct >= EPSILONEMISSIONS)
		      {
		         ballzero = false;
		      }
		   }

		   if (ballzero)
	           {
	      	      for (int ns = 0; ns < numstates; ns++)
	              {
		         emissionproducts_ni[ns] = EPSILONEMISSIONS;
		      }
		   }
		}
	     }
	  }

	  //initial probability in state s is initial probability times emission probability at first position
          double[] alpha_nt = alpha[0];
	  double dscale = 0;
	  double[] emissionproducts_nobserveindex =emissionproducts[traindataObservedIndex_nseq[0]];
 	  for (int ns = 0; ns < numstates; ns++)
          {
	      alpha_nt[ns] = probinit[ns] * emissionproducts_nobserveindex[ns];
	      dscale += alpha_nt[ns];
	  }
	  scale[0] = dscale;

	  //alpha_t(s)=P(o_0,...,o_t,x_t=s|lambda)
          //converts the alpha terms to probabilities

	  if (bscalebeta)
	  {

	     if (dscale == 0)
	     {
	        for (int ns = 0; ns < numstates; ns++)
	        {
	           //added
	           if ((alpha_nt[ns] < EPSILONSTATE)&&(emissionproducts_nobserveindex[ns]>0))//(ns < numstates-1))) 
	           {
	              alpha_nt[ns] = EPSILONSTATE;
	      	      dscale += EPSILONSTATE;
		   }
		}
	        scale[0] = dscale;
	     }

	     for (int ns = 0; ns < numstates; ns++)
	     {
                alpha_nt[ns] /= dscale;

		if ((alpha_nt[ns] < EPSILONSTATE)&&(emissionproducts_nobserveindex[ns]>0))// (ns < numstates-1))) 
	        {
	           alpha_nt[ns] = EPSILONSTATE;
		}
	     }
	  }
      	  else
	  {
             for (int ns = 0; ns < numstates; ns++)
	     {
                alpha_nt[ns] /= dscale;
	     }
	  }
	  //for (int ni = 0; ni < numstates; ni++)
          //{
          //   alpha_nt[ni] /= dscale;
	  //}
	
          //stores in coltransitionprobs the transpose of transitionprobs
          for (int ni = 0; ni < numstates; ni++)
          {
             double[] coltransitionprobs_ni = coltransitionprobs[ni];
             for (int nj = 0; nj < numstates; nj++)
	     {
	        coltransitionprobs_ni[nj] = transitionprobs[nj][ni];
	     }
	  }

          //forward step
          int numtime_nseq = numtime[nordered_nseq];
          for (int nt = 1; nt < numtime_nseq; nt++)
          {
             //the actual observed combination at position t	        
	     double[] alpha_ntm1 = alpha[nt-1];
	     alpha_nt = alpha[nt];
	      
	     dscale = 0;
	     emissionproducts_nobserveindex = emissionproducts[traindataObservedIndex_nseq[nt]];
	     for (int ns = 0; ns < numstates; ns++)
	     {
	        //going through each state		   

	        int transitionprobsnumCol_ns = transitionprobsnumCol[ns];
	        int[] transitionprobsindexCol_ns = transitionprobsindexCol[ns];
	        double[] coltransitionprobs_ns = coltransitionprobs[ns];

	        double dtempsum = 0;
                if (transitionprobsnumCol_ns < nsparsecutoff)
	        {
		    //if it is sparse enough then it is worth the extra array indirection here
	           for (int nj = 0; nj < transitionprobsnumCol_ns; nj++)
	           {
	               //for each next state computing inner sum of all previous alpha and the transition probability
	               //for all non-zero transitions into the state
			int nmappedindex = transitionprobsindexCol_ns[nj];
			dtempsum += coltransitionprobs_ns[nmappedindex]*alpha_ntm1[nmappedindex];
		   }
		}
	        else
	        {
                   for (int nj = 0; nj < numstates; nj++)
	           {
	              //for each next state computing inner sum of all previous alpha and the transition probability
	              //for all transitions into the state
		      dtempsum += coltransitionprobs_ns[nj]*alpha_ntm1[nj];
		   }
		}

                //multiply the transition sum by the emission probability
	        double dalphaval = dtempsum*emissionproducts_nobserveindex[ns];
                alpha_nt[ns] = dalphaval;
	        dscale += dalphaval;
	     }

	      //rescaling alpha
              scale[nt] = dscale;
              //scale_t(s)=P(o_0,...,o_t|lambda) summed over all states

	      if (bscalebeta)
	      {
	         if (dscale == 0)
	         {
		    for (int ns = 0; ns < numstates; ns++)
		    {
		       //added
		       if ((alpha_nt[ns] < EPSILONSTATE)&&(emissionproducts_nobserveindex[ns]>0))//(ns < numstates-1))) 
		       {
		          alpha_nt[ns] = EPSILONSTATE;
			  dscale += EPSILONSTATE;
		       }
		    }
		    scale[nt] = dscale;
		 }

	         for (int ns = 0; ns < numstates; ns++)
	         {
                    alpha_nt[ns] /= dscale;

		    if ((alpha_nt[ns] < EPSILONSTATE)&&(emissionproducts_nobserveindex[ns]>0))//(ns < numstates-1))) 
		    {
		       alpha_nt[ns] = EPSILONSTATE;
		    }
		 }
	      }
      	      else
	      {
	         for (int ns = 0; ns < numstates; ns++)
	         {
                    alpha_nt[ns] /= dscale;
		 }
	      }
	      //for (int ns = 0; ns < numstates; ns++)
              //{
	      //  alpha_nt[ns] /= dscale;
	      //}      	       
	  }
	    
          //backward step
          //beta_t(s)=P(o_t+1,...,o_T|x_t=s,lambda)
          int nlastindex = numtime_nseq-1;
	  double dinitval;
          if (bscalebeta)
	  {
	     dinitval = 1.0/numstates;
	  }
	  else
	  {
             dinitval = 1.0/scale[nlastindex];
	  }

          for (int ns = 0; ns < numstates; ns++)
	  {
	     beta_ntp1[ns] = dinitval;
	  }
          //double dinitval = 1.0/scale[nlastindex];
          //for (int ns = 0; ns < numstates; ns++)
          //{
          //    beta_ntp1[ns] = dinitval;
	  //}
	
	  int nmappedindexouter;
 
	  double ddenom = 0;	      

          //gamma_nt - P(x=S| o_0,...,o_t)
          //P(o_t+1,...,o_T|x_t=s,lambda) * P(o_0,...,o_t,x_t=s|lambda)
	  double[] gamma_nt = gamma[nlastindex]; 
          for (int ns = 0; ns < gamma_nt.length; ns++)
          {
	      double dval = alpha[nlastindex][ns]*beta_ntp1[ns];
	      ddenom += dval;
	      gamma_nt[ns] = dval;
	  }

	  if (ddenom > 0)
	  {
             for (int ns = 0; ns < gamma_nt.length; ns++)
             {
	        gamma_nt[ns] /= ddenom;
	     }
	  }


          for (int nt = nlastindex - 1; nt >= 0; nt--)
          {
	     gamma_nt = gamma[nt];
	     int ntp1 = (nt+1);
		   
	     double[] emissionproducts_ncombo_ntp1 = emissionproducts[traindataObservedIndex_nseq[ntp1]];		

	     double dsumbeta = 0;
	     double dscale_nt = scale[nt];

	     for (int ns = 0; ns < numstates; ns++)
             {
	        tempproductbetaemiss[ns] = beta_ntp1[ns]*emissionproducts_ncombo_ntp1[ns];
	     }

	     if (bscaleemissions)
	     {
	        //adding here to help numerical stability
	        double dmaxval = 0;
		for (int ns = 0; ns < numstates; ns++)
		{
	           if (tempproductbetaemiss[ns] > dmaxval)
		   {
	              dmaxval = tempproductbetaemiss[ns];
		   }
		}
                for (int ns = 0; ns < numstates; ns++)
		{
	           tempproductbetaemiss[ns] = tempproductbetaemiss[ns]/dmaxval;
		}
	     }

	     //double dscaleinv = 1.0/scale[nt];
             //scale_t(s)=P(o_0,...,o_t|lambda) summed over all states
	     for (int ni = 0; ni < numstates; ni++)
	     {
	        double dtempsum = 0;
	        int[] transitionprobsindex_ni =  transitionprobsindex[ni];
       	        double[] transitionprobs_ni = transitionprobs[ni];
		int transitionprobsnum_ni = transitionprobsnum[ni];

                if (transitionprobsnum_ni < nsparsecutoff)
	        {
	           //if it is sparse enough then it is worth the extra array indirection here
	           for (int nj = 0; nj < transitionprobsnum_ni; nj++)
	           {
	              //for each state summing over transition probability to state j, emission probablity in j at next step
	              //and probability of observing the remaining sequence
		      nmappedindexouter = transitionprobsindex_ni[nj];
		      dtempsum += transitionprobs_ni[nmappedindexouter]*tempproductbetaemiss[nmappedindexouter];			
		   }
		}
	        else
	        {
                   for (int nj = 0; nj < numstates; nj++)
	           {
	              //for each state summing over transition probability to state j, emission probablity in j at next step
	              //and probability of observing the remaining sequence
		      dtempsum += transitionprobs_ni[nj]*tempproductbetaemiss[nj];
		   }
		}

	        if (bscalebeta)
	        {
	           beta_nt[ni] = dtempsum;
       	           dsumbeta += dtempsum;
		}
	        else
		{
		   double dratio = dtempsum/dscale_nt;
		   if (dratio > Double.MAX_VALUE)
		   {
	              beta_nt[ni] = Double.MAX_VALUE;//dtempsum/dscale_nt;
		   }
		   else
		   {
		      beta_nt[ni] = dratio;
		   }
		}

		  //double dratio = dtempsum/dscale_nt;
		  //if (dratio > Double.MAX_VALUE)
		  //{
		  //    beta_nt[ni] = Double.MAX_VALUE;//dtempsum/dscale_nt;
		  //}
		  //else
		  //{
		  //    beta_nt[ni] = dratio;// dtempsum/dscale_nt;
		  //}
	      }

	      if (bscalebeta)
	      {
                 for (int ni = 0; ni < numstates; ni++)
	         {
	            beta_nt[ni]/= dsumbeta;

		    if (beta_nt[ni] < EPSILONSTATE)//&&(!bdummy))// || (ni < numstates-1))) 
		    {
		       beta_nt[ni] = EPSILONSTATE;
		    }
		 }
	      }		

	      ddenom = 0;		
	      alpha_nt = alpha[nt];

	      //gamma_nt - P(x=S| o_0,...,o_t)
              //P(o_t+1,...,o_T|x_t=s,lambda) * P(o_0,...,o_t,xt=s|lambda)

	      for (int ns = 0; ns < gamma_nt.length; ns++)
              {
	         double dval = alpha_nt[ns]*beta_nt[ns];

		 ddenom += dval;
		 gamma_nt[ns] = dval;
	      }

	      if (ddenom > 0)
	      {
	         for (int ns = 0; ns < gamma_nt.length; ns++)
                 {
	            gamma_nt[ns]/=ddenom;       		   
		 }
	      }
	      beta_ntp1 = beta_nt;		
	  }


          for (int nt = 0; nt < numtime_nseq; nt++)
	  {

             gamma_nt = gamma[nt];

	     //handling the first line
	     if ((breadsegment)||(breadstatebyline))
	     {
	        double dmaxval = 0;
                int nmaxstate = 0;

                for (int ns = 0; ns < gamma_nt.length; ns++)
                {	   	
	           double dprob = gamma_nt[ns];
	     
	           if (dprob > dmaxval)
	           {
		      //best one found so far 
	              dmaxval = dprob;
	              nmaxstate = ns;
		   }
		}

		confusion[fullmax[nt]][nmaxstate]++;
	     }
	     else
	     {
		 double[] fullposterior_nt = fullposterior[nt];
		 for (int nstate = 0; nstate < numstates; nstate++)
		 {
		     double dfullval = fullposterior_nt[nstate];
		     double dpartialval = gamma_nt[nstate];
		     if (dfullval >= dpartialval)
		     {
			 //assigned less to this state with the subset of the marks adding that amount to the decifict
			 dstatesagree[nstate] += dpartialval;
			 deficit[nstate] = dfullval - dpartialval;
			 surplus[nstate] = 0;
		     } 
		     else
		     {
			 //we have a surplus of posterior assigned to this state with a subset of marks
			 dstatesagree[nstate] += dfullval;
			 surplus[nstate] = dpartialval - dfullval;
			 deficit[nstate] = 0;
		     }
		 }

		 double dsumdenom = 0;
		 for (int nb = 0; nb < surplus.length; nb++)
		 {
		    dsumdenom += surplus[nb];
		 }
		 for (int nb = 0; nb < surplus.length; nb++)
		 {
	            //re-normalize surplus
		    surplus[nb] /= dsumdenom;
		 }

	         for (int nb = 0; nb < confusion.length; nb++)
		 {
		    double[] confusion_nb = confusion[nb];
		    if (deficit[nb] > 0)
		    {
		       double ddeficit_nb = deficit[nb];
		       for (int nc = 0; nc < confusion_nb.length; nc++)
		       {
			   confusion_nb[nc] += ddeficit_nb*surplus[nc];
                           //there is a deficit for state nb with the subset of marks
			   //allocating it to the states that proportionally have additional posterior
		       }
		    }
		 }		 		 
	     	  
	         for (int nb = 0; nb < confusion.length; nb++)
	         {
	            confusion[nb][nb] = dstatesagree[nb];
		 }
	     }
	  }
       }

       System.out.println("Writing to file "+szconfusionfileprefix+".txt");
       if (bprintimage)
       {
          System.out.println("Writing to file "+szconfusionfileprefix+".svg");
          System.out.println("Writing to file "+szconfusionfileprefix+".png");
       }
       PrintWriter pwconfusion = new PrintWriter(new FileWriter(szconfusionfileprefix+".txt",bappend));
       pwconfusion.print("EvalSubset\t"+szincludemarks);
       pwconfusion.println("\t"+szdatasets);

       for (int na = 0; na < confusion.length; na++)
       {
	   pwconfusion.print("\t"+chorder+(na+1));
       }
       pwconfusion.println();

       for (int na = 0; na < confusion.length; na++)
       {
	   pwconfusion.print(""+chorder+(na+1));
	   double ddenom = 0;
           for (int nb = 0; nb < confusion[na].length; nb++)
	   {
	       ddenom += confusion[na][nb];
	   }

	   for (int nb = 0; nb < confusion[na].length; nb++)
	   {
	       normalizedconfusion[na][nb] = confusion[na][nb]/(double) ddenom;
	       pwconfusion.print("\t"+nf.format(normalizedconfusion[na][nb]));
	   }
	   pwconfusion.println();
       }
       pwconfusion.close();
       if (bprintimage)
       {
          printConfusionImage(normalizedconfusion, szconfusionfileprefix, szincludemarks);
       }
    }


    /**
     * Takes an existing model and segmentation and outputs a confusion matrix by a selected subset
     */
    public void makeSegmentationConfusionWithLoad() throws IOException
    {
	//added in 1.15
        NumberFormat nf = NumberFormat.getInstance();
        nf.setMaximumFractionDigits(4);

       //number of non-zero transition required to be less than this at the more stringent cutoff 
       //for trying to exploit sparsity in the transition matrix for efficiency gains
       int nsparsecutoff = (int) (numstates * ChromHMM.SPARSECUTOFFRATIO);

       //int[] numtime = new int[traindataObservedIndex.length];

       //stores the maximum number of locations in any sequence and in each sequence
       int nmaxtime = 0;
       for (int nseq = 0; nseq < numtime.length; nseq++)
       {
	   //numtime[nseq] = traindataObservedIndex[nseq].length;
          if (numtime[nseq] > nmaxtime)
	  {
      	     nmaxtime = numtime[nseq];
	  }
       }

       //double
       double[][] fullposterior = null;

       int[] fullmax = null;
       double[][] confusion = new double[numstates][numstates];
       double[][] normalizedconfusion = new double[numstates][numstates];

       if (breadposterior)
       {
          fullposterior = new double[nmaxtime][numstates];
       }

       if ((breadstatebyline)||(breadsegment))
       {
	  fullmax = new int[nmaxtime];
       }


       if (ChromHMM.BVERBOSE)
       {
          System.out.println("Maximum number of locations\t"+nmaxtime);
       }


       //stores the emission probability for the i^th combination of marks in the j^th state
       double[][] emissionproducts = new double[nmaxtime][numstates]; //previously traindataObservedValues.length

       //stores temporary product terms
       double[] tempproductbetaemiss = new double[numstates];

       //This stores the alpha values at each time point and number of states
       double[][] alpha = new double[nmaxtime][numstates];

       //Temporary storage of the gamma's for each state
       double[][] gamma = new double[nmaxtime][numstates];

       //Temporary storage of the beta values for each state
       double[] beta_nt = new double[numstates];

       //Temporary storage of the beta values for each state at the next time point
       double[] beta_ntp1 = new double[numstates];

       //stores the scaling value for each time point
       double[] scale = new double[nmaxtime];

       //stores the transition probabilities for each column
       double[][] coltransitionprobs = new double[numstates][numstates];


       //stores the indicies of the data
       int[] traindataObservedIndex = new int[nmaxtime];

       //saving the mapping of signatures and chromsome observed on
       //stores whether there is a present call at each location
       boolean[][] traindataObservedValues = new boolean[nmaxtime][numdatasets];

       //stores whether the mark is not considered missing
       boolean[][] traindataNotMissing = new boolean[nmaxtime][numdatasets]; //usually nobserved

       boolean[] includemarks = new boolean[numdatasets];

       double[] surplus = new double[numstates];
       double[] deficit = new double[numstates];
       double[] dstatesagree = new double[numstates];

       if (szincludemarks.length()!=numdatasets)
       {
	   throw new IllegalArgumentException("Number of marks in "+szincludemarks+" of "+szincludemarks.length()+" does not equal expected number of "+numdatasets);
       }
      

       if ((breadstatebyline)||(breadsegment))
       {
	   //stores the maximum assignment with all marks
          fullmax = new int[nmaxtime];
       }
       else
       {
	   //stores the posterior assignment with all marks
          fullposterior = new double[nmaxtime][numstates];
       }

       String szdatasets = "";
       for (int nmark = 0; nmark < includemarks.length; nmark++)
       {
	   if (szincludemarks.charAt(nmark) == '1')
	   {
	       //stores in includemarks those data sets that have a '1' for the mark
	       includemarks[nmark] = true;
	       if (szdatasets.equals(""))
	       {
		   szdatasets += datasets[nmark];
	       }
	       else
	       {
		   szdatasets += "," + datasets[nmark];
	       }
	   }
	   else if (szincludemarks.charAt(nmark) == '0')
	   {
	       includemarks[nmark] = false;
	   }
	   else
	   {
	       throw new IllegalArgumentException(szincludemarks+" is not a valid bit string for includemarks!");
	   }
       }

       RecIntString[] ordered = new RecIntString[chromfiles.length];
       for (int nindex = 0; nindex < ordered.length; nindex++)
       {
	   ordered[nindex] = new RecIntString(nindex,chromfiles[nindex]);
       }
       Arrays.sort(ordered,new RecIntStringCompare());



       hsprefix = new HashSet();

       for (int nseq = 0; nseq <  chromfiles.length; nseq++) // traindataObservedIndex.length
       {
          int nordered_nseq = ordered[nseq].nindex;
	  //goes through each sequence

	  ///////////////////////////////////////////////////////////////////////////////
	  //load data here

	  HashMap hmObserved = new HashMap(); //maps an observation string to an index and set of flags

	  int nobserved = 0;

	  if (ChromHMM.BVERBOSE)
	  {
	     System.out.println("reading\t"+szinputdir+" "+chromfiles[nordered_nseq]);
          }
	  BufferedReader br = Util.getBufferedReader(szinputdir+"/"+chromfiles[nordered_nseq]); //fixed bug was chromfiles[nseq]
	  String szLine;
	  br.readLine(); //first line tells cell type and chromosome
	  br.readLine();//flush mark header

	  ArrayList aldata = new ArrayList();
	  while ((szLine = br.readLine())!=null)
	  {
	      StringTokenizer st = new StringTokenizer(szLine,"\t ");
	      StringBuffer sb = new StringBuffer();

	      for (int ncol = 0; ncol < numdatasets; ncol++)
	      {
	         if (!st.hasMoreTokens())
		 {
		    throw new IllegalArgumentException("Found line without "+numdatasets+" values in file "+chromfiles[nordered_nseq]);
		 }

		 String sztoken = st.nextToken();//.trim();

		 if (sztoken.equals("0"))
		 {
		    sb.append("0");
	         }
		 else if (sztoken.equals("1"))
	         {
		     sb.append("1");
		 }
		 else if (sztoken.equals("2"))
	         {
		    //this means missing
		     sb.append("2");
		 }
		 else
		 {
		    throw new IllegalArgumentException("Unrecognized value "+sztoken+" found in "+szinputdir+"/"+chromfiles[nordered_nseq]);
		 }
	      }
	      aldata.add(sb.toString());
	  }
	  br.close();

	  int nsize = aldata.size();

	  for (int nrow = 0; nrow < nsize; nrow++)
	  {
	     BigInteger theBigInteger = new BigInteger((String) aldata.get(nrow),3);
	     Integer theObservedInt  = (Integer) hmObserved.get(theBigInteger);
		  //boolean[] flagA;

	     if (theObservedInt == null)
	     {
	        //System.out.println(szmappingbyte.length());
		//storing a mapping from observed byte string to an integer index in alFlags and alObserved
		hmObserved.put(theBigInteger, Integer.valueOf(nobserved));

		//saving this observed index
		traindataObservedIndex[nrow] = nobserved;

		//increments the number of observed combinations of marks
	        nobserved++;
	     }
	     else
	     {
	        //storing the index of the flags associated with this row
		traindataObservedIndex[nrow] = ((Integer) theObservedInt).intValue();
	     }
	  }

	  Iterator hmObservedIterator = hmObserved.entrySet().iterator();
	  while (hmObservedIterator.hasNext())
	  {
	     Map.Entry pairs = (Map.Entry) hmObservedIterator.next();
	     BigInteger theBigInteger = (BigInteger) pairs.getKey();
	     String szmapping = theBigInteger.toString(3);  //getting back the mapping string

	     //ObservedRec theObservedRec = (ObservedRec) pairs.getValue();
	     int ncurrindex = ((Integer) pairs.getValue()).intValue();// theObservedRec.nobserved;//this is an index on which obervation combination it is

	     boolean[] traindataObservedValues_ncurrindex = traindataObservedValues[ncurrindex];
	     boolean[] traindataNotMissing_ncurrindex = traindataNotMissing[ncurrindex];

	     //if the mapping string is less than the number of data sets then
	     //there are leading 0's will set for leading 0's not missing and absent
	     int numch = szmapping.length();
	     int numleading0 = numdatasets - numch;
	     for (int nj = 0; nj < numleading0; nj++)
	     {
	        traindataObservedValues_ncurrindex[nj] = false;
	        traindataNotMissing_ncurrindex[nj] = true;
	     }

	     int nmappedindex = numleading0; //starting from the leading 0 position
	     for (int nj = 0; nj < numch; nj++)
	     {
		 char ch = szmapping.charAt(nj);

		 if (ch == '0')
		 {
		     traindataObservedValues_ncurrindex[nmappedindex] = false;
		     traindataNotMissing_ncurrindex[nmappedindex] = true;
		 }
		 else if (ch=='1')
		 {
		     traindataObservedValues_ncurrindex[nmappedindex] = true;
		     traindataNotMissing_ncurrindex[nmappedindex] = true;
		 }
		 else
		 {
		     //missing data
		     traindataObservedValues_ncurrindex[nmappedindex] = false;
		     traindataNotMissing_ncurrindex[nmappedindex] = false;
		 }
		 nmappedindex++;
	     }
	  }


          //int[] traindataObservedIndex_nseq = traindataObservedIndex[nordered_nseq];
          //boolean[] traindataObservedSeqFlags_nseq = traindataObservedSeqFlags[nordered_nseq];

	  String szprefix = "";
	  if (!cellSeq[nordered_nseq].equals(""))
	  {
	     szprefix += cellSeq[nordered_nseq]+"_";
	  }
	  szprefix += numstates;
	  if (!szoutfileID.equals(""))
	  {
	     szprefix += "_"+szoutfileID;
          }
	  hsprefix.add(szprefix);
	  
	  if (breadposterior)
	  {
	     BufferedReader brprobs = null;
	     //creates the posterior file
	     String szposteriorinfilename = szsegmentdir+"/POSTERIOR/"+szprefix+"_"+chromSeq[nordered_nseq]+ChromHMM.SZPOSTERIOREXTENSION;

	     //updated in v.1.11 to check for  gz extensions
	     File f = new File (szposteriorinfilename);
	     File fgz = new File (szposteriorinfilename+".gz");
	     if ((!f.exists())&&(fgz.exists()))
	     {
		 szposteriorinfilename = szposteriorinfilename + ".gz";
	     }
	     brprobs = Util.getBufferedReader(szposteriorinfilename);// new BufferedReader(new FileReader(szposteriorinfilename));

	     //skips the header lines
	     brprobs.readLine();
	     brprobs.readLine();
	     String szLinePosterior;

	     int nline = 0;
	     while ((szLinePosterior = brprobs.readLine())!=null)
	     {
		StringTokenizer stposterior = new StringTokenizer(szLinePosterior,"\t ");
	        for (int nstate = 0; nstate < numstates; nstate++)
	        {
		   fullposterior[nline][nstate] = Double.parseDouble(stposterior.nextToken());
		}
		nline++;
	     }
	     brprobs.close(); 
	  } 
          else if (breadstatebyline)
	  {
	      String szcurrchrom = chromSeq[nordered_nseq];
	      //reads a file which has the state with the maximum posterior probability
	      String szmaxinfilename = szsegmentdir+"/STATEBYLINE/"+szprefix+"_"+szcurrchrom+ChromHMM.SZSTATEBYLINEEXTENSION;

	      //updated in v.1.11 to check for  gz extensions
	      File f = new File (szmaxinfilename);
	      File fgz = new File (szmaxinfilename+".gz");
	      if ((!f.exists())&&(fgz.exists()))
	      {
	         szmaxinfilename = szmaxinfilename + ".gz";
	      }
 	      BufferedReader brmax = Util.getBufferedReader(szmaxinfilename);// new BufferedReader(new FileReader(szmaxinfilename));
	      //skip the header lines
	      brmax.readLine();
	      brmax.readLine();
	      String szLineMax;
	      int nline = 0;
	      while ((szLineMax = brmax.readLine())!=null)
	      {
		  fullmax[nline] = Integer.parseInt(szLineMax)-1;
		  nline++;
	      }
	      brmax.close(); 
	  }
	  else if (breadsegment)
	  {

	     BufferedReader brbed = null;
	     //creates a file which has the maximum segmentation
	     //we only have one file per cell type here
	     String szcurrchrom = chromSeq[nordered_nseq];

	     String szsegmentinfilename = szsegmentdir+"/" + szprefix+ChromHMM.SZSEGMENTEXTENSION;

	     //updated in v.1.11 to check for  gz extensions
	     File f = new File (szsegmentinfilename);
	     File fgz = new File (szsegmentinfilename+".gz");
	     if ((!f.exists())&&(fgz.exists()))
	     {
                szsegmentinfilename = szsegmentinfilename + ".gz";
             }

	     brbed = Util.getBufferedReader(szsegmentinfilename); // new BufferedReader(new FileReader(szsegmentinfilename));
		 
	     String szLineMax;
	     while ((szLineMax = brbed.readLine())!=null)
	     {
		 StringTokenizer stchrom = new StringTokenizer(szLineMax,"\t");
		 String szchrom = stchrom.nextToken().trim();

		 if (szchrom.equals(szcurrchrom))
	         {
		     int nbegin = Integer.parseInt(stchrom.nextToken().trim())/nbinsize;
		     int nend = (Integer.parseInt(stchrom.nextToken().trim())-1)/nbinsize;
		     int nstate = Integer.parseInt(stchrom.nextToken().trim().substring(1))-1;		        		
		     for (int nj = nbegin; nj <= nend; nj++)
		     {
			 fullmax[nj] = nstate;
		     }
		 }	      
	     }
	     brbed.close();		
	  }


          if (bscaleemissions)
	  {
	     for (int ni = 0; ni < nobserved; ni++) //emissionproducts.length
	     {
	        //going through each combination of marks
		//if (traindataObservedSeqFlags_nseq[ni])
		//{
	        //this signature of marks is observed on the current chromosome so
		//updating its emission probabilities
                double[] emissionproducts_ni = emissionproducts[ni];
		boolean[] traindataObservedValues_ni = traindataObservedValues[ni];
		boolean[] traindataNotMissing_ni = traindataNotMissing[ni];

		for (int ns = 0; ns < numstates; ns++)
		{
		   emissionproducts_ni[ns] = 1;
		}

		for (int nmod = 0; nmod < numdatasets; nmod++)
		{
		   for (int ns = 0; ns < numstates; ns++)
		   {
		      if ((traindataNotMissing_ni[nmod])&&(includemarks[nmod]))
		      {
		         //we are include this marks emission probability
			 if (traindataObservedValues_ni[nmod])
			 {
		            emissionproducts_ni[ns] *= emissionprobs[ns][nmod][1];
			 }
		         else
		         {
		            emissionproducts_ni[ns] *= emissionprobs[ns][nmod][0];
			 }
		      }
			 // otherwise treated as missing omitting from product
		   }

		   double dmaxval = 0;
		   for (int ns = 0; ns < numstates; ns++)
		   {
		      if (emissionproducts_ni[ns] > dmaxval)
		      {
		         dmaxval = emissionproducts_ni[ns];
		      }
		   }

		   if (dmaxval <= 0)
		   {
		      for (int ns = 0; ns < numstates; ns++)
		      {
		         emissionproducts_ni[ns] = 1;
		      }
		   }
		   else
		   {
		      for (int ns = 0; ns < numstates; ns++)
		      {
		         emissionproducts_ni[ns]/= dmaxval;
		      }
		   }
		}
	     }
	     //}
	  }
	  else
	  {
	     for (int ni = 0; ni < nobserved; ni++)
             {
	        //going through each combination of marks
	        //if (traindataObservedSeqFlags_nseq[ni])
	        //{
	        //this signature of marks is observed on the current chromosome so
	        //updating its emission probabilities
	        double[] emissionproducts_ni = emissionproducts[ni];
	        boolean[] traindataObservedValues_ni = traindataObservedValues[ni];
	        boolean[] traindataNotMissing_ni = traindataNotMissing[ni];	
	        boolean ballzero = true;
	  
	        for (int ns = 0; ns < numstates; ns++)
	        {
	           double dproduct = 1;
	           double[][] emissionprobs_ni = emissionprobs[ns];

		   //going through all marks
		   for (int nmod = 0; nmod < numdatasets; nmod++)
	           {
		      if ((traindataNotMissing_ni[nmod])&&(includemarks[nmod]))
		      {
		         //we have observed the mark
		         if (traindataObservedValues_ni[nmod])
		         {
		            dproduct *= emissionprobs_ni[nmod][1];
		         }
		         else 
	                 {
		            dproduct *= emissionprobs_ni[nmod][0];
		         }
		      }
		         // otherwise treated as missing omitting from product
		   }
	           emissionproducts_ni[ns] = dproduct;

		   if (dproduct >= EPSILONEMISSIONS)
		   {
	              ballzero = false;
		   }
		}

		if (ballzero)
	        {
	      	   for (int ns = 0; ns < numstates; ns++)
	           {
		      emissionproducts_ni[ns] = EPSILONEMISSIONS;
		   }
		}
	     }
	  }

	  //initial probability in state s is initial probability times emission probability at first position
          double[] alpha_nt = alpha[0];
	  double dscale = 0;
	  double[] emissionproducts_nobserveindex =emissionproducts[traindataObservedIndex[0]];
 	  for (int ns = 0; ns < numstates; ns++)
          {
	      alpha_nt[ns] = probinit[ns] * emissionproducts_nobserveindex[ns];
	      dscale += alpha_nt[ns];
	  }
	  scale[0] = dscale;

	  //alpha_t(s)=P(o_0,...,o_t,x_t=s|lambda)
          //converts the alpha terms to probabilities
	  if (bscalebeta)
          {
	     if (dscale == 0)
	     {
	        for (int ns = 0; ns < numstates; ns++)
		{
	           //added
		   if ((alpha_nt[ns] < EPSILONSTATE)&&(emissionproducts_nobserveindex[ns]>0))//(ns < numstates-1))) 
	           {
	              alpha_nt[ns] = EPSILONSTATE;
		      dscale += EPSILONSTATE;
		   }
		}
	        scale[0] = dscale;
	     }

             for (int ns = 0; ns < numstates; ns++)
	     {
                alpha_nt[ns] /= dscale;

	        if ((alpha_nt[ns] < EPSILONSTATE)&& (emissionproducts_nobserveindex[ns]>0))//(ns < numstates-1))) 
	        {
	           alpha_nt[ns] = EPSILONSTATE;
	        }
             }
	  }
      	  else
	  {
	     for (int ns = 0; ns < numstates; ns++)
	     {
                alpha_nt[ns] /= dscale;
	     }
	  }
	  //for (int ni = 0; ni < numstates; ni++)
          //{
          //   alpha_nt[ni] /= dscale;
	  //}
	
          //stores in coltransitionprobs the transpose of transitionprobs
          for (int ni = 0; ni < numstates; ni++)
          {
             double[] coltransitionprobs_ni = coltransitionprobs[ni];
             for (int nj = 0; nj < numstates; nj++)
	     {
	        coltransitionprobs_ni[nj] = transitionprobs[nj][ni];
	     }
	  }

          //forward step
          int numtime_nseq = numtime[nordered_nseq];
          for (int nt = 1; nt < numtime_nseq; nt++)
          {
             //the actual observed combination at position t	        
	     double[] alpha_ntm1 = alpha[nt-1];
	     alpha_nt = alpha[nt];
	      
	     dscale = 0;
	     emissionproducts_nobserveindex = emissionproducts[traindataObservedIndex[nt]];
	     for (int ns = 0; ns < numstates; ns++)
	     {
	        //going through each state		   

	        int transitionprobsnumCol_ns = transitionprobsnumCol[ns];
	        int[] transitionprobsindexCol_ns = transitionprobsindexCol[ns];
	        double[] coltransitionprobs_ns = coltransitionprobs[ns];

	        double dtempsum = 0;
                if (transitionprobsnumCol_ns < nsparsecutoff)
	        {
		    //if it is sparse enough then it is worth the extra array indirection here
	           for (int nj = 0; nj < transitionprobsnumCol_ns; nj++)
	           {
	               //for each next state computing inner sum of all previous alpha and the transition probability
	               //for all non-zero transitions into the state
			int nmappedindex = transitionprobsindexCol_ns[nj];
			dtempsum += coltransitionprobs_ns[nmappedindex]*alpha_ntm1[nmappedindex];
		   }
		}
	        else
	        {
                   for (int nj = 0; nj < numstates; nj++)
	           {
	              //for each next state computing inner sum of all previous alpha and the transition probability
	              //for all transitions into the state
		      dtempsum += coltransitionprobs_ns[nj]*alpha_ntm1[nj];
		   }
		}

                //multiply the transition sum by the emission probability
	        double dalphaval = dtempsum*emissionproducts_nobserveindex[ns];
                alpha_nt[ns] = dalphaval;
	        dscale += dalphaval;
	     }

	      //rescaling alpha
              scale[nt] = dscale;
              //scale_t(s)=P(o_0,...,o_t|lambda) summed over all states

	      if (bscalebeta)
	      {
	         if (dscale == 0)
	         {
		    for (int ns = 0; ns < numstates; ns++)
		    {
		       //added
		       if ((alpha_nt[ns] < EPSILONSTATE)&&(emissionproducts_nobserveindex[ns]>0))//(ns < numstates-1))) 
		       {
		          alpha_nt[ns] = EPSILONSTATE;
			  dscale += EPSILONSTATE;
		       }
		    }
		    scale[nt] = dscale;
		 }

	         for (int ns = 0; ns < numstates; ns++)
	         {
                    alpha_nt[ns] /= dscale;

		    if ((alpha_nt[ns] < EPSILONSTATE)&&(emissionproducts_nobserveindex[ns]>0))//(ns < numstates-1))) 
		    {
		       alpha_nt[ns] = EPSILONSTATE;
		    }
		 }
	      }
      	      else
	      {
	         for (int ns = 0; ns < numstates; ns++)
	         {
                    alpha_nt[ns] /= dscale;
		 }
	      }
	      //for (int ns = 0; ns < numstates; ns++)
              //{
	      //  alpha_nt[ns] /= dscale;
	      //}      	       
	  }
	    
          //backward step
          //beta_t(s)=P(o_t+1,...,o_T|x_t=s,lambda)
          int nlastindex = numtime_nseq-1;
	  double dinitval;
	  if (bscalebeta)
	  {
	     dinitval = 1.0/numstates;
	  }
	  else
	  {
             dinitval = 1.0/scale[nlastindex];
	  }

          for (int ns = 0; ns < numstates; ns++)
	  {
             beta_ntp1[ns] = dinitval;
	  }
          //double dinitval = 1.0/scale[nlastindex];
          //for (int ns = 0; ns < numstates; ns++)
          //{
          //    beta_ntp1[ns] = dinitval;
	  //}
	
	  int nmappedindexouter;
 
	  double ddenom = 0;	      

          //gamma_nt - P(x=S| o_0,...,o_t)
          //P(o_t+1,...,o_T|x_t=s,lambda) * P(o_0,...,o_t,x_t=s|lambda)
	  double[] gamma_nt = gamma[nlastindex]; 
          for (int ns = 0; ns < gamma_nt.length; ns++)
          {
	      double dval = alpha[nlastindex][ns]*beta_ntp1[ns];
	      ddenom += dval;
	      gamma_nt[ns] = dval;
	  }

	  if (ddenom > 0)
	  {
             for (int ns = 0; ns < gamma_nt.length; ns++)
             {
	        gamma_nt[ns] /= ddenom;
	     }
	  }


          for (int nt = nlastindex - 1; nt >= 0; nt--)
          {
	      gamma_nt = gamma[nt];
	      int ntp1 = (nt+1);
		   
	      double[] emissionproducts_ncombo_ntp1 = emissionproducts[traindataObservedIndex[ntp1]];		
	      double dsumbeta = 0;
	      double dscale_nt = scale[nt];

	      for (int ns = 0; ns < numstates; ns++)
              {
		  tempproductbetaemiss[ns] = beta_ntp1[ns]*emissionproducts_ncombo_ntp1[ns];
	      }

	      if (bscaleemissions)
	      {
		  //adding here to help numerical stability
		  double dmaxval = 0;
		  for (int ns = 0; ns < numstates; ns++)
		  {
		      if (tempproductbetaemiss[ns] > dmaxval)
		      {
			  dmaxval = tempproductbetaemiss[ns];
		      }
		  }
                  for (int ns = 0; ns < numstates; ns++)
		  {
		     tempproductbetaemiss[ns] = tempproductbetaemiss[ns]/dmaxval;
		  }
	      }

	      //double dscaleinv = 1.0/scale[nt];
              //scale_t(s)=P(o_0,...,o_t|lambda) summed over all states
	      for (int ni = 0; ni < numstates; ni++)
	      {
		  double dtempsum = 0;
		  int[] transitionprobsindex_ni =  transitionprobsindex[ni];
		  double[] transitionprobs_ni = transitionprobs[ni];
		  int transitionprobsnum_ni = transitionprobsnum[ni];

                  if (transitionprobsnum_ni < nsparsecutoff)
	          {
		    //if it is sparse enough then it is worth the extra array indirection here
	             for (int nj = 0; nj < transitionprobsnum_ni; nj++)
	             {
	                //for each state summing over transition probability to state j, emission probablity in j at next step
	                //and probability of observing the remaining sequence
		        nmappedindexouter = transitionprobsindex_ni[nj];
		        dtempsum += transitionprobs_ni[nmappedindexouter]*tempproductbetaemiss[nmappedindexouter];			
		     }
		  }
	          else
	          {
                     for (int nj = 0; nj < numstates; nj++)
	             {
	                //for each state summing over transition probability to state j, emission probablity in j at next step
	                //and probability of observing the remaining sequence
		        dtempsum += transitionprobs_ni[nj]*tempproductbetaemiss[nj];
		     }
		  }

		  if (bscalebeta)
		  {
		     beta_nt[ni] = dtempsum;
		     dsumbeta += dtempsum;
		  }
		  else
		  {
		     double dratio = dtempsum/dscale_nt;
		     if (dratio > Double.MAX_VALUE)
		     {
		        beta_nt[ni] = Double.MAX_VALUE;//dtempsum/dscale_nt;
		     }
		     else
		     {
		        beta_nt[ni] = dratio;
		     }
		  }
		  //double dratio = dtempsum/dscale_nt;
		  //if (dratio > Double.MAX_VALUE)
		  //{
		  //    beta_nt[ni] = Double.MAX_VALUE;//dtempsum/dscale_nt;
		  //}
		  //else
		  //{
		  //    beta_nt[ni] = dratio;// dtempsum/dscale_nt;
		  //}
	      }

	      if (bscalebeta)
	      {
                 for (int ni = 0; ni < numstates; ni++)
	         {
	            beta_nt[ni]/= dsumbeta;

		    if (beta_nt[ni] < EPSILONSTATE)//&&(!bdummy))// || (ni < numstates-1))) 
		    {
		       beta_nt[ni] = EPSILONSTATE;
		    }
		 }
	      }		

	      ddenom = 0;		
	      alpha_nt = alpha[nt];

	       //gamma_nt - P(x=S| o_0,...,o_t)
               //P(o_t+1,...,o_T|x_t=s,lambda) * P(o_0,...,o_t,xt=s|lambda)

	       for (int ns = 0; ns < gamma_nt.length; ns++)
               {
		   double dval = alpha_nt[ns]*beta_nt[ns];

		   ddenom += dval;
		   gamma_nt[ns] = dval;
	       }

	       if (ddenom > 0)
	       {
	          for (int ns = 0; ns < gamma_nt.length; ns++)
                  {
		     gamma_nt[ns]/=ddenom;       		   
		  }
	       }
	       beta_ntp1 = beta_nt;		
	  }


          for (int nt = 0; nt < numtime_nseq; nt++)
	  {

             gamma_nt = gamma[nt];

	     //handling the first line
	     if ((breadsegment)||(breadstatebyline))
	     {
	        double dmaxval = 0;
                int nmaxstate = 0;

                for (int ns = 0; ns < gamma_nt.length; ns++)
                {	   	
	           double dprob = gamma_nt[ns];
	     
	           if (dprob > dmaxval)
	           {
		      //best one found so far 
	              dmaxval = dprob;
	              nmaxstate = ns;
		   }
		}

		confusion[fullmax[nt]][nmaxstate]++;
	     }
	     else
	     {
		 double[] fullposterior_nt = fullposterior[nt];
		 for (int nstate = 0; nstate < numstates; nstate++)
		 {
		     double dfullval = fullposterior_nt[nstate];
		     double dpartialval = gamma_nt[nstate];
		     if (dfullval >= dpartialval)
		     {
			 //assigned less to this state with the subset of the marks adding that amount to the decifict
			 dstatesagree[nstate] += dpartialval;
			 deficit[nstate] = dfullval - dpartialval;
			 surplus[nstate] = 0;
		     } 
		     else
		     {
			 //we have a surplus of posterior assigned to this state with a subset of marks
			 dstatesagree[nstate] += dfullval;
			 surplus[nstate] = dpartialval - dfullval;
			 deficit[nstate] = 0;
		     }
		 }

		 double dsumdenom = 0;
		 for (int nb = 0; nb < surplus.length; nb++)
		 {
		    dsumdenom += surplus[nb];
		 }
		 for (int nb = 0; nb < surplus.length; nb++)
		 {
	            //re-normalize surplus
		    surplus[nb] /= dsumdenom;
		 }

	         for (int nb = 0; nb < confusion.length; nb++)
		 {
		    double[] confusion_nb = confusion[nb];
		    if (deficit[nb] > 0)
		    {
		       double ddeficit_nb = deficit[nb];
		       for (int nc = 0; nc < confusion_nb.length; nc++)
		       {
			   confusion_nb[nc] += ddeficit_nb*surplus[nc];
                           //there is a deficit for state nb with the subset of marks
			   //allocating it to the states that proportionally have additional posterior
		       }
		    }
		 }		 		 
	     	  
	         for (int nb = 0; nb < confusion.length; nb++)
	         {
	            confusion[nb][nb] = dstatesagree[nb];
		 }
	     }
	  }
       }

       System.out.println("Writing to file "+szconfusionfileprefix+".txt");
       if (bprintimage)
       {
          System.out.println("Writing to file "+szconfusionfileprefix+".svg");
          System.out.println("Writing to file "+szconfusionfileprefix+".png");
       }
       PrintWriter pwconfusion = new PrintWriter(new FileWriter(szconfusionfileprefix+".txt",bappend));
       pwconfusion.print("EvalSubset\t"+szincludemarks);
       pwconfusion.println("\t"+szdatasets);

       for (int na = 0; na < confusion.length; na++)
       {
	   pwconfusion.print("\t"+chorder+(na+1));
       }
       pwconfusion.println();

       for (int na = 0; na < confusion.length; na++)
       {
	   pwconfusion.print(""+chorder+(na+1));
	   double ddenom = 0;
           for (int nb = 0; nb < confusion[na].length; nb++)
	   {
	       ddenom += confusion[na][nb];
	   }

	   for (int nb = 0; nb < confusion[na].length; nb++)
	   {
	       normalizedconfusion[na][nb] = confusion[na][nb]/(double) ddenom;
	       pwconfusion.print("\t"+nf.format(normalizedconfusion[na][nb]));
	   }
	   pwconfusion.println();
       }
       pwconfusion.close();
       if (bprintimage)
       {
          printConfusionImage(normalizedconfusion, szconfusionfileprefix, szincludemarks);
       }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////

    /**
     * Outputs the confusion matrix as a '.png' and '.svg' 
     */
    public void printConfusionImage(double[][] confusion, String szconfusionfileprefix,
				    String szincludemarks) throws IOException
    {

        String[] rowlabels = new String[numstates];

        for (int ni = 0; ni < numstates; ni++)
        {
	    rowlabels[ni] = ""+(ni+1);
	    String szsuffix;

	    //changing in v1.12 to be state labels for new state to match char ordering
	    //if ((szsuffix = (String) hmlabelExtend.get(""+chorder+(stateordering[ni]+1)))!=null)
	    if ((szsuffix = (String) hmlabelExtend.get(""+chorder+(ni+1)))!=null)
	    {
		rowlabels[ni]+= "_"+szsuffix;
	    }
	}


        HeatChart map = new HeatChart(confusion);

        map.setTitle("Confusion Matrix");
        map.setXAxisLabel("State Subset of Marks ("+szincludemarks+")");
        map.setAxisValuesFont(new Font("SansSerif",0,20));
        map.setAxisLabelsFont(new Font("SansSerif",0,22));
        map.setTitleFont(new Font("SansSerif",0,24));
        map.setYAxisLabel("State All Marks");
        if (confusion.length <=5)
        {
	   map.setChartMargin(125);
        }
        else
        {
	   map.setChartMargin(100);
        }
        map.setXValues(rowlabels);
        map.setYValues(rowlabels);
        map.setHighValueColour(theColor);


        Util.printImageToSVG(map, szconfusionfileprefix+".svg");
        map.saveToFile(new File(szconfusionfileprefix+".png"));

    }


    //////////////////////////////////////////////////////////////////////////////////////////

    /**
     * Takes an existing model and outputs information about the segmentation depending on the values of
     * bprintsegment, bprintstatebyline, bprintposterior
     */
    public void makeSegmentationWithLoad() throws IOException
    {
        NumberFormat nf = NumberFormat.getInstance(Locale.ENGLISH); //updated in 1.14 to Local.ENGLISH to ensure format of posterior values can be read in 
        nf.setMaximumFractionDigits(4);

       //number of non-zero transition required to be less than this at the more stringent cutoff 
       //for trying to exploit sparsity in the transition matrix for efficiency gains
       int nsparsecutoff = (int) (numstates * ChromHMM.SPARSECUTOFFRATIO);

       //int[] numtime = new int[traindataObservedIndex.length];

       //stores the maximum number of locations in any sequence and in each sequence
       int nmaxtime = 0;
       for (int nseq = 0; nseq < numtime.length; nseq++)
       {
	   //numtime[nseq] = traindataObservedIndex[nseq].length;
          if (numtime[nseq] > nmaxtime)
	  {
      	     nmaxtime = numtime[nseq];
	  }
       }


       HashMap hmMaxCoord = null;
       if (szchromlengthfile != null)
       {
	   hmMaxCoord = new HashMap();
	   BufferedReader brchromlengthfile =  Util.getBufferedReader(szchromlengthfile);
	   String szLine;
	   while ((szLine = brchromlengthfile.readLine())!=null)
	   {
	       StringTokenizer st = new StringTokenizer(szLine,"\t ");
	       hmMaxCoord.put(st.nextToken(),Integer.valueOf(st.nextToken()));
	   }
	   brchromlengthfile.close();
       }

       if (ChromHMM.BVERBOSE)
       {
          System.out.println("Maximum number of locations\t"+nmaxtime);
       }

       //stores the emission probability for the i^th combination of marks in the j^th state
       double[][] emissionproducts = new double[nmaxtime][numstates]; //previously traindataObservedValues.length

       //stores temporary product terms
       double[] tempproductbetaemiss = new double[numstates];

       //This stores the alpha values at each time point and number of states
       double[][] alpha = new double[nmaxtime][numstates];

       //Temporary storage of the gamma's for each state
       double[][] gamma = new double[nmaxtime][numstates];

       //Temporary storage of the beta values for each state
       double[] beta_nt = new double[numstates];

       //Temporary storage of the beta values for each state at the next time point
       double[] beta_ntp1 = new double[numstates];

       //stores the scaling value for each time point
       double[] scale = new double[nmaxtime];

       //stores the transition probabilities for each column
       double[][] coltransitionprobs = new double[numstates][numstates];
       
       //stores the indicies of the data
       int[] traindataObservedIndex = new int[nmaxtime];

       //saving the mapping of signatures and chromsome observed on
       //stores whether there is a present call at each location
       boolean[][] traindataObservedValues = new boolean[nmaxtime][numdatasets];

       //stores whether the mark is not considered missing
       boolean[][] traindataNotMissing = new boolean[nmaxtime][numdatasets]; //usually nobserved


       //maps cell ID to printwriter objects
       HashMap hmcellToFourColPW = null;
       if (bprintsegment)
       {
          hmcellToFourColPW = new HashMap();
       }

       RecIntString[] ordered = null;//new RecIntString[chromfiles.length];
       RecIntStringSplit[] orderedsplit = null;//new RecIntString[chromfiles.length];

       if (bsplit)
       {
          orderedsplit = new RecIntStringSplit[chromfiles.length];
	  for (int nindex = 0; nindex < orderedsplit.length; nindex++)
          {
	     int nlastperiodindex = chromSeq[nindex].lastIndexOf('.');

	     if (nlastperiodindex == -1)
	     {
	        throw new IllegalArgumentException("No period found in chromosome "+chromSeq[nindex]+" despite split being specified");
	     }

	     String szchromportion = chromSeq[nindex].substring(0, nlastperiodindex);
	     int nsplitbinindex = Integer.parseInt(chromSeq[nindex].substring(nlastperiodindex+1));
             orderedsplit[nindex] = new RecIntStringSplit(nindex,cellSeq[nindex],szchromportion,nsplitbinindex);
	  }
	  Arrays.sort(orderedsplit,new RecIntStringSplitCompare());
       }
       else
       {
          ordered = new RecIntString[chromfiles.length];
	  for (int nindex = 0; nindex < ordered.length; nindex++)
          {
	     ordered[nindex] = new RecIntString(nindex,chromfiles[nindex]);
	  }
	  Arrays.sort(ordered,new RecIntStringCompare());
       }


       //RecIntString[] ordered = new RecIntString[chromfiles.length];
       //for (int nindex = 0; nindex < ordered.length; nindex++)
       //{
       //	   ordered[nindex] = new RecIntString(nindex,chromfiles[nindex]);
       //}
       //Arrays.sort(ordered,new RecIntStringCompare());

       hsprefix = new HashSet();

       GZIPOutputStream pwprobszip = null;
       GZIPOutputStream pwmaxzip = null;
       GZIPOutputStream pwbedzip = null;
       PrintWriter pwprobs = null;
       PrintWriter pwmax = null;
       PrintWriter pwbed = null;

       int noffset = 0;
       boolean bclosefile = true;
       boolean bnewfile = true;


       //these can go across split files so defined out here
       int nstart = 0; //the start index of the current active interval
       int nmaxstateprev = -1;
       int nprevlinecount = 0;
       for (int nseq = 0; nseq < chromfiles.length; nseq++)
       {
	   int nordered_nseq;// = ordered[nseq].nindex;
           String szactualchrom;
           if (bsplit)
	   {
	      nordered_nseq = orderedsplit[nseq].nindex;
	      szactualchrom = orderedsplit[nseq].szchrom;
           }
           else
           {
	      nordered_nseq = ordered[nseq].nindex;
	      szactualchrom = chromSeq[nordered_nseq];
	   }

	   if (bsplit)
	   {
	       //noffset = orderedsplit[nseq].nsplitbinindex * numsplitbins;
	       //if (noffset == 0)
	      if (orderedsplit[nseq].nsplitbinindex == 0)
	      {
		 noffset = 0;
	         bnewfile = true;
	      }
	      else
	      {
		 noffset += nprevlinecount;
	         if ((nseq >= 1) && (orderedsplit[nseq].nsplitbinindex != (orderedsplit[nseq-1].nsplitbinindex+1)))
		 {
		    throw new IllegalArgumentException("For "+orderedsplit[nseq].szcell+"_"+orderedsplit[nseq].szchrom+" found "+
                                                         "a file with split index "+orderedsplit[nseq].nsplitbinindex+", but not "+
								      (orderedsplit[nseq].nsplitbinindex-1));
		 }
	         bnewfile = false;
	      }
	      bclosefile = ((nseq + 1 == orderedsplit.length)||(orderedsplit[nseq+1].nsplitbinindex == 0));
	   }


	   //goes through each sequence

	   ///////////////////////////////////////////////////////////////////////////////
	   //load data here

	   HashMap hmObserved = new HashMap(); //maps an observation string to an index and set of flags
 
    	   int nobserved = 0;

	   if (ChromHMM.BVERBOSE)
	   {
              System.out.println("reading\t"+szinputdir+" "+chromfiles[nordered_nseq]);
	   }
	   BufferedReader br = Util.getBufferedReader(szinputdir+"/"+chromfiles[nordered_nseq]); //fixed bug was chromfiles[nseq]
           String szLine; 
           br.readLine(); //first line tells cell type and chromosome
	   br.readLine();//flush mark header
	   ArrayList aldata = new ArrayList();
	   while ((szLine = br.readLine())!=null)
	   {
	      StringTokenizer st = new StringTokenizer(szLine,"\t ");
      	      StringBuffer sb = new StringBuffer();
		
	      for (int ncol = 0; ncol < numdatasets; ncol++)
	      {
	         if (!st.hasMoreTokens())
		 {
	            throw new IllegalArgumentException("Found line without "+numdatasets+" values in file "+chromfiles[nordered_nseq]);
		 }

		 String sztoken = st.nextToken();
		    
		 if (sztoken.equals("0"))
	         {
	            sb.append("0");
	         }
		 else if (sztoken.equals("1"))
	         {
	            sb.append("1");
		 }
	         else if (sztoken.equals("2"))
		 {
	            //this means missing
	            sb.append("2");
		 }
		 else
	         {
	            throw new IllegalArgumentException("Unrecognized value "+sztoken+" found in "+szinputdir+"/"+chromfiles[nordered_nseq]);
	         }
	      }
	      aldata.add(sb.toString());
	   }
	   br.close();
	      
	   int nsize = aldata.size();
	   nprevlinecount = nsize;
	   //traindataObservedIndex[nfile] = new int[nsize];
	   //int[] traindataObservedIndex_nfile = traindataObservedIndex[nfile];

	   for (int nrow = 0; nrow < nsize; nrow++)
	   {
	      BigInteger theBigInteger = new BigInteger((String) aldata.get(nrow),3);
	      Integer theObservedInt  = (Integer) hmObserved.get(theBigInteger);
	       //boolean[] flagA;

	      if (theObservedInt == null)
       	      {
	         //System.out.println(szmappingbyte.length());
	         //storing a mapping from observed byte string to an integer index in alFlags and alObserved
	         hmObserved.put(theBigInteger, Integer.valueOf(nobserved));

		 //saving this observed index
		 traindataObservedIndex[nrow] = nobserved;

		 //increments the number of observed combinations of marks
		 nobserved++;
	      }
	      else
	      {
	         //storing the index of the flags associated with this row 
      	         traindataObservedIndex[nrow] = ((Integer) theObservedInt).intValue();
	      }
	   }  	 
	   
	   Iterator hmObservedIterator = hmObserved.entrySet().iterator();
	   while (hmObservedIterator.hasNext())
	   {
	      Map.Entry pairs = (Map.Entry) hmObservedIterator.next();
	      BigInteger theBigInteger = (BigInteger) pairs.getKey();
	      String szmapping = theBigInteger.toString(3);  //getting back the mapping string

	      //ObservedRec theObservedRec = (ObservedRec) pairs.getValue();
	      int ncurrindex = ((Integer) pairs.getValue()).intValue();// theObservedRec.nobserved;//this is an index on which obervation combination it is

	      boolean[] traindataObservedValues_ncurrindex = traindataObservedValues[ncurrindex];
	      boolean[] traindataNotMissing_ncurrindex = traindataNotMissing[ncurrindex]; 
	   
	      //if the mapping string is less than the number of data sets then 
	      //there are leading 0's will set for leading 0's not missing and absent
	      int numch = szmapping.length();
	      int numleading0 = numdatasets - numch;
	      for (int nj = 0; nj < numleading0; nj++)
	      {
	         traindataObservedValues_ncurrindex[nj] = false;
	         traindataNotMissing_ncurrindex[nj] = true;
	      }

	      int nmappedindex = numleading0; //starting from the leading 0 position
	      for (int nj = 0; nj < numch; nj++)
	      {
	         char ch = szmapping.charAt(nj);

	         if (ch == '0')
	         {
		    traindataObservedValues_ncurrindex[nmappedindex] = false;
		    traindataNotMissing_ncurrindex[nmappedindex] = true;
		 }
	         else if (ch=='1')
	         {
		    traindataObservedValues_ncurrindex[nmappedindex] = true;
		    traindataNotMissing_ncurrindex[nmappedindex] = true;
		 }
	         else
	         {
		      //missing data
		    traindataObservedValues_ncurrindex[nmappedindex] = false;
		    traindataNotMissing_ncurrindex[nmappedindex] = false;
		 }
	         nmappedindex++;
	      }
	   }
       
	   //////////////////////////////////////////////////////////////////////////////

	   //int[] traindataObservedIndex_nseq = traindataObservedIndex[nordered_nseq];
           //boolean[] traindataObservedSeqFlags_nseq = traindataObservedSeqFlags[nordered_nseq];

	   if (bnewfile)
	   {
	      String szprefix = "";
	      if (!cellSeq[nordered_nseq].equals(""))
	      {
	         szprefix += cellSeq[nordered_nseq]+"_";
	      }
	      szprefix += numstates;
	      if (!szoutfileID.equals(""))
	      {
	         szprefix += "_"+szoutfileID;
	      }
	      hsprefix.add(szprefix);

	      //GZIPOutputStream pwprobszip = null;
	      //PrintWriter pwprobs = null;

	      //GZIPOutputStream pwmaxzip = null;
	      //PrintWriter pwmax = null;

              //GZIPOutputStream pwbedzip = null;
	      //PrintWriter pwbed = null;

	      if (bgzip)
	      {
	         //PrintWriter pwprobs = null;
	         if (bprintposterior)
	         {
	            //creates the posterior file
 
		     //changed in v.1.18
	            //String szposterioroutfilename = szoutputdir+"/POSTERIOR/"+szprefix+"_"+chromSeq[nordered_nseq]+ChromHMM.SZPOSTERIOREXTENSION+".gz";
	            String szposterioroutfilename = szoutputdir+"/POSTERIOR/"+szprefix+"_"+szactualchrom+ChromHMM.SZPOSTERIOREXTENSION+".gz";

                    System.out.println("Writing to file "+szposterioroutfilename);
		    pwprobszip = new GZIPOutputStream(new FileOutputStream(szposterioroutfilename));
	            //pwprobs = new PrintWriter(szposterioroutfilename);
		    String szout = cellSeq[nordered_nseq]+"\t"+szactualchrom+"\n";

		    byte[] btformat = szout.getBytes();
		    pwprobszip.write(btformat,0,btformat.length);

		    StringBuffer sbout = new StringBuffer();
	            //pwprobs.println(cellSeq[nordered_nseq]+"\t"+chromSeq[nordered_nseq]);
	            for (int ni = 0; ni < numstates-1; ni++)
	            {
		       sbout.append(""+chorder+(ni+1)+"\t");
		       //pwprobs.print(""+chorder+(ni+1)+"\t");
		    } 
		    sbout.append(""+chorder+(numstates)+"\n");
	            //pwprobs.println(""+chorder+(numstates));
		    btformat = sbout.toString().getBytes();
		    pwprobszip.write(btformat,0,btformat.length);
		 }

   	         //PrintWriter pwmax = null;

	         if (bprintstatebyline)
	         {
	            //creates a file which has the state with the maximum posterior probability
	            String szmaxoutfilename = szoutputdir+"/STATEBYLINE/"+szprefix+"_"+szactualchrom+ChromHMM.SZSTATEBYLINEEXTENSION+".gz";

	            System.out.println("Writing to file "+szmaxoutfilename);
	            //pwmax = new PrintWriter(szmaxoutfilename);
		    pwmaxzip = new GZIPOutputStream(new FileOutputStream(szmaxoutfilename));

		    String szout = cellSeq[nordered_nseq]+"\t"+szactualchrom+"\n";
                    byte[] btformat = szout.getBytes();
                    pwmaxzip.write(btformat,0,btformat.length); //fix 1.23
	            //pwmax.println(cellSeq[nordered_nseq]+"\t"+chromSeq[nordered_nseq]);

		    szout = "MaxState "+chorder+"\n";
		    btformat = szout.getBytes();
                    pwmaxzip.write(btformat,0,btformat.length); //fix 1.23
	            //pwmax.println("MaxState "+chorder);
		 }

	         //PrintWriter pwbed = null;
	         //GZIPOutputStream pwbedzip = null;

	         if (bprintsegment)
	         {
	            //creates a file which has the maximum segmentation
	            //we only have one file per cell type here
	            pwbedzip = (GZIPOutputStream) hmcellToFourColPW.get(cellSeq[nordered_nseq]);
		    // (PrintWriter) hmcellToFourColPW.get(cellSeq[nordered_nseq]);

	            if (pwbedzip == null)
	            {
		       //haven't seen this cell type
	               String szsegmentoutfilename = szoutputdir+"/" + szprefix+SZSEGMENTEXTENSION+".gz";

		       pwbedzip = new GZIPOutputStream(new FileOutputStream(szsegmentoutfilename));
	               //pwbed = new PrintWriter(szsegmentoutfilename);
		       System.out.println("Writing to file "+szsegmentoutfilename);

	               hmcellToFourColPW.put(cellSeq[nordered_nseq],pwbedzip);    
		    }
		 }
	      }
	      else
	      {

	         if (bprintposterior)
	         {
	            //creates the posterior file
 
	            String szposterioroutfilename = szoutputdir+"/POSTERIOR/"+szprefix+"_"+szactualchrom+ChromHMM.SZPOSTERIOREXTENSION;

                    System.out.println("Writing to file "+szposterioroutfilename);
	            pwprobs = new PrintWriter(szposterioroutfilename);
	            pwprobs.println(cellSeq[nordered_nseq]+"\t"+szactualchrom);
	            for (int ni = 0; ni < numstates-1; ni++)
	            {
		       pwprobs.print(""+chorder+(ni+1)+"\t");
		    }
	            pwprobs.println(""+chorder+(numstates));
		 }

	         if (bprintstatebyline)
	         {
	            //creates a file which has the state with the maximum posterior probability
	            String szmaxoutfilename = szoutputdir+"/STATEBYLINE/"+szprefix+"_"+szactualchrom+ChromHMM.SZSTATEBYLINEEXTENSION;

	            System.out.println("Writing to file "+szmaxoutfilename);
	            pwmax = new PrintWriter(szmaxoutfilename);
	            pwmax.println(cellSeq[nordered_nseq]+"\t"+szactualchrom);
	            pwmax.println("MaxState "+chorder);
		 }

	         //PrintWriter pwbed = null;
	         if (bprintsegment)
	         {
	            //creates a file which has the maximum segmentation
	            //we only have one file per cell type here
	            pwbed = (PrintWriter) hmcellToFourColPW.get(cellSeq[nordered_nseq]);

	            if (pwbed == null)
	            {
		       //haven't seen this cell type
	               String szsegmentoutfilename = szoutputdir+"/" + szprefix+SZSEGMENTEXTENSION;

	               pwbed = new PrintWriter(szsegmentoutfilename);
		       System.out.println("Writing to file "+szsegmentoutfilename);

	               hmcellToFourColPW.put(cellSeq[nordered_nseq],pwbed);    
		    }
		 }
	      }
	   }


           if (bscaleemissions)
  	   {
	      for (int ni = 0; ni < nobserved; ni++)
	      {
	         //going through each combination of marks
	         //if (traindataObservedSeqFlags_nseq[ni])
	         //{
	 	 //this signature of marks is observed on the current chromosome so
		 //updating its emission probabilities
		 double[] emissionproducts_ni = emissionproducts[ni];
	         boolean[] traindataObservedValues_ni = traindataObservedValues[ni];
		 boolean[] traindataNotMissing_ni = traindataNotMissing[ni];

		 for (int ns = 0; ns < numstates; ns++)
		 {
		    emissionproducts_ni[ns] = 1;
		 }

		 for (int nmod = 0; nmod < numdatasets; nmod++)
		 {
		    for (int ns = 0; ns < numstates; ns++)
		    {
		       if (traindataNotMissing_ni[nmod])
		       {
		          //we are include this marks emission probability
			  if (traindataObservedValues_ni[nmod])
			  {
			     emissionproducts_ni[ns] *= emissionprobs[ns][nmod][1];
			  }
			  else
			  {
			     emissionproducts_ni[ns] *= emissionprobs[ns][nmod][0];
			  }
                       }
		       // otherwise treated as missing omitting from product
		    }

                    double dmaxval = 0;
                    for (int ns = 0; ns < numstates; ns++)
		    {
		       if (emissionproducts_ni[ns] > dmaxval)
		       {
		          dmaxval = emissionproducts_ni[ns];
		       }
		    }


                    if (dmaxval <= 0)//EPSILONEMISSIONS)
		    {
		       for (int ns = 0; ns < numstates; ns++)
		       {
		          emissionproducts_ni[ns] = 1;
		       }
		    }
		    else
		    {
		       for (int ns = 0; ns < numstates; ns++)
		       {
		          emissionproducts_ni[ns]/= dmaxval;
		       }
		    }
		 }
	      }
	   }
	   else
	   {
	      for (int ni = 0; ni < nobserved; ni++)
              {
	         //going through each combination of marks
	         //if (traindataObservedSeqFlags_nseq[ni])
	         //{
	         //this signature of marks is observed on the current chromosome so
	         //updating its emission probabilities
	         double[] emissionproducts_ni = emissionproducts[ni];
	         boolean[] traindataObservedValues_ni = traindataObservedValues[ni];
	         boolean[] traindataNotMissing_ni = traindataNotMissing[ni];		  

	         boolean ballzero = true;

	         for (int ns = 0; ns < numstates; ns++)
	         {
	            double dproduct = 1;
	            double[][] emissionprobs_ni = emissionprobs[ns];

		    //going through all marks
	            for (int nmod = 0; nmod < numdatasets; nmod++)
	            {
		       if (traindataNotMissing_ni[nmod])
		       {
		          //we have observed the mark
		          if (traindataObservedValues_ni[nmod])
		          {
		             dproduct *= emissionprobs_ni[nmod][1];
		          }
		          else 
	                  {
		             dproduct *= emissionprobs_ni[nmod][0];
		          }
		       }
		       // otherwise treated as missing omitting from product
		    }
	            emissionproducts_ni[ns] = dproduct;

		    if (dproduct >= EPSILONEMISSIONS)
		    {
	               ballzero = false;
		    }
		 }
	     
	         if (ballzero)
	         {
	            for (int ns = 0; ns < numstates; ns++)
	            {
	               emissionproducts_ni[ns] = EPSILONEMISSIONS;
		    }
		 }
	      }
	   }
	  
	   //initial probability in state s is initial probability times emission probability at first position
           double[] alpha_nt = alpha[0];
	   double dscale = 0;
	   double[] emissionproducts_nobserveindex =emissionproducts[traindataObservedIndex[0]];
 	   for (int ns = 0; ns < numstates; ns++)
           {
	      alpha_nt[ns] = probinit[ns] * emissionproducts_nobserveindex[ns];
	      dscale += alpha_nt[ns];
	   }
	   scale[0] = dscale;

	   //alpha_t(s)=P(o_0,...,o_t,x_t=s|lambda)
           //converts the alpha terms to probabilities
	   if (bscalebeta)
	   {
	      if (dscale == 0)
	      {
	         for (int ns = 0; ns < numstates; ns++)
		 {
	             //added
		     if ((alpha_nt[ns] < EPSILONSTATE)&&(emissionproducts_nobserveindex[ns]>0))//(ns < numstates-1))) 
		     {
		        alpha_nt[ns] = EPSILONSTATE;
		        dscale += EPSILONSTATE;
		     }
		 }
	         scale[0] = dscale;
	      }

	      for (int ns = 0; ns < numstates; ns++)
	      {
                 alpha_nt[ns] /= dscale;

	         if ((alpha_nt[ns] < EPSILONSTATE)&&(emissionproducts_nobserveindex[ns]>0))//(ns < numstates-1))) 
       	         {
	            alpha_nt[ns] = EPSILONSTATE;
		 }
	      }
	   }
      	   else
           {
	      for (int ns = 0; ns < numstates; ns++)
	      {
                 alpha_nt[ns] /= dscale;
	      }
	   }
	   //for (int ni = 0; ni < numstates; ni++)
           //{
           //   alpha_nt[ni] /= dscale;
	   //}
	
           //stores in coltransitionprobs the transpose of transitionprobs
           for (int ni = 0; ni < numstates; ni++)
           {
              double[] coltransitionprobs_ni = coltransitionprobs[ni];
              for (int nj = 0; nj < numstates; nj++)
	      {
	         coltransitionprobs_ni[nj] = transitionprobs[nj][ni];
	      }
	   }

           //forward step
           int numtime_nseq = numtime[nordered_nseq];
           for (int nt = 1; nt < numtime_nseq; nt++)
           {
              //the actual observed combination at position t	        
	      double[] alpha_ntm1 = alpha[nt-1];
	      alpha_nt = alpha[nt];
	      
	      dscale = 0;
	      emissionproducts_nobserveindex = emissionproducts[traindataObservedIndex[nt]];
	      for (int ns = 0; ns < numstates; ns++)
	      {
	         //going through each state		   
	         int transitionprobsnumCol_ns = transitionprobsnumCol[ns];
	         int[] transitionprobsindexCol_ns = transitionprobsindexCol[ns];
	         double[] coltransitionprobs_ns = coltransitionprobs[ns];

	         double dtempsum = 0;
                 if (transitionprobsnumCol_ns < nsparsecutoff)
	         {
		    //if it is sparse enough then it is worth the extra array indirection here
	            for (int nj = 0; nj < transitionprobsnumCol_ns; nj++)
	            {
	               //for each next state computing inner sum of all previous alpha and the transition probability
	               //for all non-zero transitions into the state
	               int nmappedindex = transitionprobsindexCol_ns[nj];
	       	       dtempsum += coltransitionprobs_ns[nmappedindex]*alpha_ntm1[nmappedindex];
		    }
		 }
	         else
	         {
                    for (int nj = 0; nj < numstates; nj++)
	            {
	               //for each next state computing inner sum of all previous alpha and the transition probability
	               //for all transitions into the state
		       dtempsum += coltransitionprobs_ns[nj]*alpha_ntm1[nj];
		    }
		 }

                 //multiply the transition sum by the emission probability
	         double dalphaval = dtempsum*emissionproducts_nobserveindex[ns];
                 alpha_nt[ns] = dalphaval;
	         dscale += dalphaval;
	      }

	      //rescaling alpha
              scale[nt] = dscale;
              //scale_t(s)=P(o_0,...,o_t|lambda) summed over all states

	      if (bscalebeta)
	      {
	         if (dscale == 0)
	         {
		    for (int ns = 0; ns < numstates; ns++)
		    {
		       //added
		       if ((alpha_nt[ns] < EPSILONSTATE)&&(emissionproducts_nobserveindex[ns]>0))//(ns < numstates-1))) 
		       {
		          alpha_nt[ns] = EPSILONSTATE;
			  dscale += EPSILONSTATE;
		       }
		    }
		    scale[nt] = dscale;
		 }

	         for (int ns = 0; ns < numstates; ns++)
	         {
                    alpha_nt[ns] /= dscale;

		    if ((alpha_nt[ns] < EPSILONSTATE)&&(emissionproducts_nobserveindex[ns]>0))//(ns < numstates-1))) 
		    {
		       alpha_nt[ns] = EPSILONSTATE;
		    }
		 }
	      }
      	      else
	      {
	         for (int ns = 0; ns < numstates; ns++)
	         {
                    alpha_nt[ns] /= dscale;
		 }
	      }
	     //for (int ns = 0; ns < numstates; ns++)
	     //{
	     //	  alpha_nt[ns] /= dscale;
	     //}      	       
	  }
	    
          //backward step
          //beta_t(s)=P(o_t+1,...,o_T|x_t=s,lambda)
          int nlastindex = numtime_nseq-1;
	  double dinitval;
	  if (bscalebeta)
	  {
	     dinitval = 1.0/numstates;
	  }
	  else
	  {
             dinitval = 1.0/scale[nlastindex];
	  }

          for (int ns = 0; ns < numstates; ns++)
          {
             beta_ntp1[ns] = dinitval;
          }
          //double dinitval = 1.0/scale[nlastindex];
          //for (int ns = 0; ns < numstates; ns++)
          //{
          //    beta_ntp1[ns] = dinitval;
	  //}
	
	  int nmappedindexouter;
 
	  double ddenom = 0;	      

          //gamma_nt - P(x=S| o_0,...,o_t)
          //P(o_t+1,...,o_T|x_t=s,lambda) * P(o_0,...,o_t,x_t=s|lambda)
	  double[] gamma_nt = gamma[nlastindex]; 
          for (int ns = 0; ns < gamma_nt.length; ns++)
          {
	      double dval = alpha[nlastindex][ns]*beta_ntp1[ns];
	      ddenom += dval;
	      gamma_nt[ns] = dval;
	  }

	  if (ddenom > 0)
	  {
             for (int ns = 0; ns < gamma_nt.length; ns++)
             {
	        gamma_nt[ns] /= ddenom;
	     }
	  }


          for (int nt = nlastindex - 1; nt >= 0; nt--)
          {
	      gamma_nt = gamma[nt];
	      int ntp1 = (nt+1);
		   
	      double[] emissionproducts_ncombo_ntp1 = emissionproducts[traindataObservedIndex[ntp1]];		
	      double dsumbeta = 0;
	      double dscale_nt = scale[nt];

	      for (int ns = 0; ns < numstates; ns++)
              {
		  tempproductbetaemiss[ns] = beta_ntp1[ns]*emissionproducts_ncombo_ntp1[ns];
	      }

	      if (bscaleemissions)
	      {
		  //adding here to help numerical stability
		  double dmaxval = 0;
		  for (int ns = 0; ns < numstates; ns++)
		  {
		      if (tempproductbetaemiss[ns] > dmaxval)
		      {
			  dmaxval = tempproductbetaemiss[ns];
		      }
		  }
                  for (int ns = 0; ns < numstates; ns++)
		  {
		      tempproductbetaemiss[ns] = tempproductbetaemiss[ns]/dmaxval;
		  }
	      }

	      //double dscaleinv = 1.0/scale[nt];
              //scale_t(s)=P(o_0,...,o_t|lambda) summed over all states
	      for (int ni = 0; ni < numstates; ni++)
	      {
		  double dtempsum = 0;
		  int[] transitionprobsindex_ni =  transitionprobsindex[ni];
		  double[] transitionprobs_ni = transitionprobs[ni];
		  int transitionprobsnum_ni = transitionprobsnum[ni];

                  if (transitionprobsnum_ni < nsparsecutoff)
	          {
		    //if it is sparse enough then it is worth the extra array indirection here
	             for (int nj = 0; nj < transitionprobsnum_ni; nj++)
	             {
	                //for each state summing over transition probability to state j, emission probablity in j at next step
	                //and probability of observing the remaining sequence
		        nmappedindexouter = transitionprobsindex_ni[nj];
		        dtempsum += transitionprobs_ni[nmappedindexouter]*tempproductbetaemiss[nmappedindexouter];			
		     }
		  }
	          else
	          {
                     for (int nj = 0; nj < numstates; nj++)
	             {
	                //for each state summing over transition probability to state j, emission probablity in j at next step
	                //and probability of observing the remaining sequence
		        dtempsum += transitionprobs_ni[nj]*tempproductbetaemiss[nj];
		     }
		  }

		  if (bscalebeta)
		  {
		     beta_nt[ni] = dtempsum;
		     dsumbeta += dtempsum;
		  }
		  else
		  {
		     double dratio = dtempsum/dscale_nt;
		     if (dratio > Double.MAX_VALUE)
		     {
		        beta_nt[ni] = Double.MAX_VALUE;//dtempsum/dscale_nt;
		     }
		     else
		     {
		        beta_nt[ni] = dratio;
		     }
		  }
		   //double dratio = dtempsum/dscale_nt;
		   //if (dratio > Double.MAX_VALUE)
		   //{
		   //   beta_nt[ni] = Double.MAX_VALUE;
		   //}
		  //else
		  //{
		  //  beta_nt[ni] = dratio;//dtempsum/dscale_nt;
		  //}
	      }

	      if (bscalebeta)
	      {
                 for (int ni = 0; ni < numstates; ni++)
	         {
	            beta_nt[ni]/= dsumbeta;

		    if (beta_nt[ni] < EPSILONSTATE)//&&(!bdummy))// || (emissionproducts_nobserveindex[ns]>0)))//(ni < numstates-1))) 
		    {
		       beta_nt[ni] = EPSILONSTATE;
		    }
		 }
	      }		

	      ddenom = 0;		
	      alpha_nt = alpha[nt];

	       //gamma_nt - P(x=S| o_0,...,o_t)
               //P(o_t+1,...,o_T|x_t=s,lambda) * P(o_0,...,o_t,xt=s|lambda)

	       for (int ns = 0; ns < gamma_nt.length; ns++)
               {
		   double dval = alpha_nt[ns]*beta_nt[ns];

		   ddenom += dval;
		   gamma_nt[ns] = dval;
	       }

	       if (ddenom > 0)
	       {
	          for (int ns = 0; ns < gamma_nt.length; ns++)
                  {
		      gamma_nt[ns]/=ddenom;       		   
		  }
	       }
	       beta_ntp1 = beta_nt;		
	  }

	  //int nstart = 0; //the start index of the current active interval
          gamma_nt = gamma[0];

	  double dmaxval = 0;
          int nmaxstate = 0;

	  int nt = 0;

	  if (bgzip)
	  {

	     if (bnewfile)
	     {
	        //handling the first line
                for (int ns = 0; ns < gamma_nt.length; ns++)
                {	   	
	           int nmappedstate = stateordering[ns]; //maps new state to old
	           double dprob = gamma_nt[nmappedstate];
	           if (bprintposterior)
	           {
		      String szout;
	              if (ns > 0)
	              {
		         //print with tab if not the first
		         szout = "\t"+nf.format(dprob);
                         //pwprobs.print("\t"+nf.format(dprob));
		      }
                      else
                      {
		         szout = nf.format(dprob);
		         //pwprobs.print(nf.format(dprob));
		      }
		      byte[] btformat = szout.getBytes();
		      pwprobszip.write(btformat,0,btformat.length);
		   }

	           if (dprob > dmaxval)
	           {
		      //best one found so far 
	              dmaxval = dprob;
	              nmaxstate = ns;
		   }
		}

	        if (bprintposterior)
	        {
		   String szout = "\n";
		   byte[] btformat = szout.getBytes();
		   pwprobszip.write(btformat,0,btformat.length);
		   //pwprobs.println();
		}

	        if (bprintstatebyline)
	        {
                   String szout = (nmaxstate+1)+"\n";
                   byte[] btformat = szout.getBytes();
                   pwmaxzip.write(btformat,0,btformat.length);
		   //pwmax.println(""+(nmaxstate+1));
	        }

                nmaxstateprev = nmaxstate;
                nt = 1;
                nstart = 0;
	     }

	     //this contains the best state of the previous interval
	     //int nmaxstateprev = nmaxstate;  	

             for (; nt < numtime_nseq; nt++)
             {
                gamma_nt = gamma[nt];

	        dmaxval = 0;
	        nmaxstate = 0;
                for (int ns = 0; ns < gamma_nt.length; ns++)
                {	   	
		   int nmappedstate = stateordering[ns]; //maps new state to old
		   double dprob = gamma_nt[nmappedstate];
 	           if (bprintposterior)
	           {
		      String szout;
	              if (ns > 0)
	              {
		         //print with tab the first time
			 szout = "\t"+nf.format(dprob);
                         //pwprobs.print("\t"+nf.format(dprob));
		      }
	              else
	              {
			 szout = nf.format(dprob);
                         //pwprobs.print(nf.format(dprob));
		      }
		      byte[] btformat = szout.getBytes();
		      pwprobszip.write(btformat,0,btformat.length);
		   }

	           if (dprob > dmaxval)
	           {
	              dmaxval = dprob;
	              nmaxstate = ns;
		   }
	        }

	        if (bprintposterior)
	        {
		   String szout = "\n";
		   byte[] btformat = szout.getBytes();
		   pwprobszip.write(btformat,0,btformat.length);
		    //pwprobs.println();
	        }

                if (bprintstatebyline)
	        {
		   String szout =""+(nmaxstate+1)+"\n";
		   byte[] btformat = szout.getBytes();
		   pwmaxzip.write(btformat,0,btformat.length);
		   //pwmax.println(""+(nmaxstate+1));  	
	        }	     
	     
	        if (bprintsegment&&(nmaxstateprev != nmaxstate))
	        {
		   //print out last segment we are done with
		   //pwbed.println(chromSeq[nordered_nseq]+"\t"+(nstart*nbinsize)+"\t"+(nt*nbinsize)+"\t"+chorder+(nmaxstateprev+1));		 
		   String szout = szactualchrom+"\t"+(nstart*nbinsize)+"\t"+((nt+noffset)*nbinsize)+"\t"+chorder+(nmaxstateprev+1)+"\n";
                   byte[] btformat = szout.getBytes();
                   pwbedzip.write(btformat,0,btformat.length);
		   //start a new segment now
		   nstart = nt+noffset;
		   nmaxstateprev = nmaxstate;
		}	     
	     }

	     if (bclosefile)
	     {
	        if (bprintsegment)
	        {
	           int nlastcoordinate;
	           Integer objMaxCoord = null;
	           if (hmMaxCoord != null)
	           {
		      objMaxCoord = ((Integer) hmMaxCoord.get(szactualchrom));
		   }

	           if (objMaxCoord != null)
	           {
		       nlastcoordinate = Math.min((numtime_nseq+noffset)*nbinsize,((Integer) objMaxCoord).intValue());

	  	      if (nstart*nbinsize >= nlastcoordinate)
		      {
			  //v1.26
			  throw new IllegalArgumentException("Start of interval "+(nstart*nbinsize+1)+" (1-based) is past the end of the chromosome "+nlastcoordinate+
                                                             "; Check if chromosome length file given to LearnModel is the same used for the binarization");
		      }
	           }
	           else
	           {
		       nlastcoordinate = (numtime_nseq+noffset)*nbinsize;
		   }

		   String szout = szactualchrom+"\t"+(nstart*nbinsize)+"\t"+nlastcoordinate+"\t"+chorder+(nmaxstateprev+1)+"\n";
		   byte[] btformat = szout.getBytes();
		   pwbedzip.write(btformat,0,btformat.length);

	        //pwbed.println(chromSeq[nordered_nseq]+"\t"+(nstart*nbinsize)+"\t"+nlastcoordinate+"\t"+chorder+(nmaxstateprev+1));
		}

	        //close out the max state file if that was requested
	        if (bprintstatebyline)
	        {
		   pwmaxzip.finish();
	           pwmaxzip.close();
	        }
	        //close out the posterior state file if that was requested
	        if (bprintposterior)
	        {
		   pwprobszip.finish();
	           pwprobszip.close();    
		}
	     }
	  }
	  else
	  {
	     if (bnewfile)
	     {
	        //handling the first line
                for (int ns = 0; ns < gamma_nt.length; ns++)
                {	   	
	           int nmappedstate = stateordering[ns]; //maps new state to old
	           double dprob = gamma_nt[nmappedstate];

	           if (bprintposterior)
	           {
	              if (ns > 0)
	              {
                         pwprobs.print("\t"+nf.format(dprob));
	              }
                      else
                      {
                         pwprobs.print(nf.format(dprob));
		      }
		   }

	           if (dprob > dmaxval)
	           {
		      //best one found so far 
	              dmaxval = dprob;
	              nmaxstate = ns;
		   }
		}

	        if (bprintposterior)
	        {
                   pwprobs.println();
	        }

	        if (bprintstatebyline)
	        {
	           pwmax.println(""+(nmaxstate+1));
	        }

	        //this contains the best state of the previous interval
	        nmaxstateprev = nmaxstate;  	
                nt = 1;
                nstart = 0;
	     }

             for ( ; nt < numtime_nseq; nt++)
             {
                gamma_nt = gamma[nt];

	        dmaxval = 0;
	        nmaxstate = 0;
                for (int ns = 0; ns < gamma_nt.length; ns++)
                {	   	
		   int nmappedstate = stateordering[ns]; //maps new state to old
		   double dprob = gamma_nt[nmappedstate];
 	           if (bprintposterior)
	           {
	              if (ns > 0)
	              {
		         //print with tab the first time
                         pwprobs.print("\t"+nf.format(dprob));
		      }
	              else
	              {
                         pwprobs.print(nf.format(dprob));
		      }
		   }

	           if (dprob > dmaxval)
	           {
	              dmaxval = dprob;
	              nmaxstate = ns;
		   }
	        }

	        if (bprintposterior)
	        {
	           pwprobs.println();
	        }

                if (bprintstatebyline)
	        {
		   pwmax.println(""+(nmaxstate+1));  	
	        }	     
	     
	        if (bprintsegment&&(nmaxstateprev != nmaxstate))
	        {
		   //print out last segment we are done with
		    pwbed.println(szactualchrom+"\t"+(nstart*nbinsize)+"\t"+((nt+noffset)*nbinsize)+"\t"+chorder+(nmaxstateprev+1));		 
		   //start a new segment now
		   nstart = nt+noffset;
		   nmaxstateprev = nmaxstate;
		}	     
	     }

	     if (bclosefile)
	     {
	        if (bprintsegment)
	        {
	           int nlastcoordinate;
	           Integer objMaxCoord = null;
	           if (hmMaxCoord != null)
	           {
		      objMaxCoord = ((Integer) hmMaxCoord.get(szactualchrom));
	           }

	           if (objMaxCoord != null)
	           {
		       nlastcoordinate = Math.min((numtime_nseq+noffset)*nbinsize,((Integer) objMaxCoord).intValue());

	  	      if (nstart*nbinsize >= nlastcoordinate)
		      {
			  //v1.26
			  throw new IllegalArgumentException("Start of interval "+(nstart*nbinsize+1)+" (1-based) is past the end of the chromosome "+nlastcoordinate+
                                                             "; Check if chromosome length file given to LearnModel is the same used for the binarization");
		      }
	           }
	           else
	           {
		       nlastcoordinate = (numtime_nseq+noffset)*nbinsize;
	           }
	           pwbed.println(szactualchrom+"\t"+(nstart*nbinsize)+"\t"+nlastcoordinate+"\t"+chorder+(nmaxstateprev+1));
		}

	        //close out the max state file if that was requested
	        if (bprintstatebyline)
	        {
	           pwmax.close();
	        }
	        //close out the posterior state file if that was requested
	        if (bprintposterior)
	        {
	           pwprobs.close();    
	        }
	     }
	  }
       }
	
       //if segment print was requested then we are going to go close those printwriters
       if (bprintsegment)
       {
          Iterator itr =  hmcellToFourColPW.values().iterator();
	  if (bgzip)
	  {
             while (itr.hasNext())
             {
	        GZIPOutputStream pwzip = (GZIPOutputStream) itr.next();
	        //PrintWriter pw = (PrintWriter) itr.next();
	        pwzip.finish();
	        pwzip.close();
	           //pw.close();
	     }
	  }
	  else
	  {
             while (itr.hasNext())
             {
                PrintWriter pw = (PrintWriter) itr.next();
	        pw.close();
	     }
	  }
       }	        
    }



    ///////////////////////////////////////////////////////////////////////////////////////////

    /**
     * Takes an existing model and outputs information about the segmentation depending on the values of
     * bprintsegment, bprintstatebyline, bprintposterior
     */
    public void makeSegmentation() throws IOException
    {
        NumberFormat nf = NumberFormat.getInstance(Locale.ENGLISH); //updated in 1.14 to Local.ENGLISH to ensure format of posterior values can be read in
        nf.setMaximumFractionDigits(4);

       //number of non-zero transition required to be less than this at the more stringent cutoff 
       //for trying to exploit sparsity in the transition matrix for efficiency gains
       int nsparsecutoff = (int) (numstates * ChromHMM.SPARSECUTOFFRATIO);

       int[] numtime = new int[traindataObservedIndex.length];

       //stores the maximum number of locations in any sequence and in each sequence
       int nmaxtime = 0;
       for (int nseq = 0; nseq < traindataObservedIndex.length; nseq++)
       {
          numtime[nseq] = traindataObservedIndex[nseq].length;
          if (numtime[nseq] > nmaxtime)
	  {
      	     nmaxtime = numtime[nseq];
	  }
       }


       HashMap hmMaxCoord = null;
       if (szchromlengthfile != null)
       {
	   hmMaxCoord = new HashMap();
	   BufferedReader brchromlengthfile =  Util.getBufferedReader(szchromlengthfile);
	   String szLine;
	   while ((szLine = brchromlengthfile.readLine())!=null)
	   {
	       StringTokenizer st = new StringTokenizer(szLine,"\t ");
	       hmMaxCoord.put(st.nextToken(),Integer.valueOf(st.nextToken()));
	   }
	   brchromlengthfile.close();
       }

       if (ChromHMM.BVERBOSE)
       {
          System.out.println("Maximum number of locations\t"+nmaxtime);
       }

       //stores the emission probability for the i^th combination of marks in the j^th state
       double[][] emissionproducts = new double[traindataObservedValues.length][numstates];

       //stores temporary product terms
       double[] tempproductbetaemiss = new double[numstates];

       //This stores the alpha values at each time point and number of states
       double[][] alpha = new double[nmaxtime][numstates];

       //Temporary storage of the gamma's for each state
       double[][] gamma = new double[nmaxtime][numstates];

       //Temporary storage of the beta values for each state
       double[] beta_nt = new double[numstates];

       //Temporary storage of the beta values for each state at the next time point
       double[] beta_ntp1 = new double[numstates];

       //stores the scaling value for each time point
       double[] scale = new double[nmaxtime];

       //stores the transition probabilities for each column
       double[][] coltransitionprobs = new double[numstates][numstates];
       
       //maps cell ID to printwriter objects
       HashMap hmcellToFourColPW = null;
       if (bprintsegment)
       {
          hmcellToFourColPW = new HashMap();
       }

       RecIntString[] ordered = null;//new RecIntString[chromfiles.length];
       RecIntStringSplit[] orderedsplit = null;//new RecIntString[chromfiles.length];

       if (bsplit)
       {
	   orderedsplit = new RecIntStringSplit[chromfiles.length];
	   for (int nindex = 0; nindex < orderedsplit.length; nindex++)
	   {
	      int nlastperiodindex = chromSeq[nindex].lastIndexOf('.');

	      if (nlastperiodindex == -1)
	      {
	         throw new IllegalArgumentException("No period found in chromosome "+chromSeq[nindex]+" despite split being specified");
	      }

	      String szchromportion = chromSeq[nindex].substring(0, nlastperiodindex);
	      int nsplitbinindex = Integer.parseInt(chromSeq[nindex].substring(nlastperiodindex+1));
	      orderedsplit[nindex] = new RecIntStringSplit(nindex,cellSeq[nindex],szchromportion,nsplitbinindex);
	   }
	   Arrays.sort(orderedsplit,new RecIntStringSplitCompare());
       }
       else
       {
          ordered = new RecIntString[chromfiles.length];
          for (int nindex = 0; nindex < ordered.length; nindex++)
          {
	     ordered[nindex] = new RecIntString(nindex,chromfiles[nindex]);
          }
          Arrays.sort(ordered,new RecIntStringCompare());
       }

       hsprefix = new HashSet();

       GZIPOutputStream pwprobszip = null;
       GZIPOutputStream pwmaxzip = null;
       GZIPOutputStream pwbedzip = null;
       PrintWriter pwprobs = null;
       PrintWriter pwmax = null;
       PrintWriter pwbed = null;

       int noffset = 0;
       boolean bclosefile = true;
       boolean bnewfile = true;

       //these can go across split files so defined out here
       int nstart = 0; //the start index of the current active interval
       int nmaxstateprev = -1;
       int nprevlinecount = 0;
       for (int nseq = 0; nseq < traindataObservedIndex.length; nseq++)
       {
	   int nordered_nseq;
	   String szactualchrom;
	   if (bsplit)
	   {
              nordered_nseq = orderedsplit[nseq].nindex;
	      szactualchrom = orderedsplit[nseq].szchrom;
	   }
	   else
	   {
	      nordered_nseq = ordered[nseq].nindex;
	      szactualchrom = chromSeq[nordered_nseq];
	   }
	   //goes through each sequence

          int[] traindataObservedIndex_nseq = traindataObservedIndex[nordered_nseq];
          boolean[] traindataObservedSeqFlags_nseq = traindataObservedSeqFlags[nordered_nseq];


	  if (bsplit)
	  {
	      //noffset = orderedsplit[nseq].nsplitbinindex * numsplitbins;
	      if (orderedsplit[nseq].nsplitbinindex==0) //(noffset == 0)
	      {
		  noffset = 0;
		  bnewfile = true; 
	      }
	      else
	      {
		  noffset += nprevlinecount;
		  if ((nseq >= 1) && (orderedsplit[nseq].nsplitbinindex != (orderedsplit[nseq-1].nsplitbinindex+1)))
		  {
		      throw new IllegalArgumentException("For "+orderedsplit[nseq].szcell+"_"+orderedsplit[nseq].szchrom+" found "+
                                                         "a file with split index "+orderedsplit[nseq].nsplitbinindex+", but not "+
                                                         (orderedsplit[nseq].nsplitbinindex-1));
		  }
		  bnewfile = false;
	      }
	      bclosefile = ((nseq + 1 == orderedsplit.length)||(orderedsplit[nseq+1].nsplitbinindex == 0));
	  }


	  if (bnewfile)
	  {

	     String szprefix = "";

	     if (!cellSeq[nordered_nseq].equals(""))
	     {
	        szprefix += cellSeq[nordered_nseq]+"_";
	     }
	     szprefix += numstates;
	     if (!szoutfileID.equals(""))
	     {
	        szprefix += "_"+szoutfileID;
             }
	     hsprefix.add(szprefix);


	     if (bgzip)
	     {
	        if (bprintposterior)
	        {
	           //creates the posterior file
		    //v.1.18 replacing chromSeq[nordered_nseq] with szactualchrom
 	           //String szposterioroutfilename = szoutputdir+"/POSTERIOR/"+szprefix+"_"+chromSeq[nordered_nseq]+ChromHMM.SZPOSTERIOREXTENSION+".gz";
	           String szposterioroutfilename = szoutputdir+"/POSTERIOR/"+szprefix+"_"+szactualchrom+ChromHMM.SZPOSTERIOREXTENSION+".gz";

                   System.out.println("Writing to file "+szposterioroutfilename);
	           //pwprobs = new PrintWriter(szposterioroutfilename);
		   pwprobszip = new GZIPOutputStream(new FileOutputStream(szposterioroutfilename));
		   String szout = cellSeq[nordered_nseq]+"\t"+szactualchrom+"\n";

		   byte[] btformat = szout.getBytes();
		   pwprobszip.write(btformat,0,btformat.length);

	           //pwprobs.println(cellSeq[nordered_nseq]+"\t"+chromSeq[nordered_nseq]);
		   StringBuffer sbout = new StringBuffer();
	           for (int ni = 0; ni < numstates-1; ni++)
	           {
		      sbout.append(""+chorder+(ni+1)+"\t");
		      //pwprobs.print(""+chorder+(ni+1)+"\t");
	           } 
		   sbout.append(""+chorder+(numstates)+"\n");

		   btformat = sbout.toString().getBytes();
                   pwprobszip.write(btformat,0,btformat.length);
	           //pwprobs.println(""+chorder+(numstates));
		}

	        if (bprintstatebyline)
	        {
	           //creates a file which has the state with the maximum posterior probability
	           String szmaxoutfilename = szoutputdir+"/STATEBYLINE/"+szprefix+"_"+szactualchrom+ChromHMM.SZSTATEBYLINEEXTENSION+".gz";
	           System.out.println("Writing to file "+szmaxoutfilename);
	 	   pwmaxzip = new GZIPOutputStream(new FileOutputStream(szmaxoutfilename));
	           //pwmax = new PrintWriter(szmaxoutfilename);
		   String szout = cellSeq[nordered_nseq]+"\t"+szactualchrom +"\n";
                   byte[] btformat = szout.getBytes();
                   pwmaxzip.write(btformat,0,btformat.length);

		   szout = "MaxState "+chorder+"\n";

	 	   btformat = szout.getBytes();
                   pwmaxzip.write(btformat,0,btformat.length);
	           //pwmax.println(cellSeq[nordered_nseq]+"\t"+chromSeq[nordered_nseq]);
	           //pwmax.println("MaxState "+chorder);
		}

	        if (bprintsegment)
	        {
	           //creates a file which has the maximum segmentation
	           //we only have one file per cell type here
	           pwbedzip = (GZIPOutputStream) hmcellToFourColPW.get(cellSeq[nordered_nseq]);

	           if (pwbedzip == null)
	           {
		      //haven't seen this cell type
	              String szsegmentoutfilename = szoutputdir+"/" + szprefix+SZSEGMENTEXTENSION+".gz";

		      pwbedzip = new GZIPOutputStream(new FileOutputStream(szsegmentoutfilename));
		      //pwbed = new PrintWriter(szsegmentoutfilename);
		      System.out.println("Writing to file "+szsegmentoutfilename);
	              hmcellToFourColPW.put(cellSeq[nordered_nseq],pwbedzip);    
		   }
		}
	     }
	     else
	     {
	        if (bprintposterior)
	        {
	           //creates the posterior file
 
	           String szposterioroutfilename = szoutputdir+"/POSTERIOR/"+szprefix+"_"+szactualchrom+ChromHMM.SZPOSTERIOREXTENSION;

                   System.out.println("Writing to file "+szposterioroutfilename);
	           pwprobs = new PrintWriter(szposterioroutfilename);
	           pwprobs.println(cellSeq[nordered_nseq]+"\t"+chromSeq[nordered_nseq]);
	           for (int ni = 0; ni < numstates-1; ni++)
	           {
		      pwprobs.print(""+chorder+(ni+1)+"\t");
	           } 
	           pwprobs.println(""+chorder+(numstates));
		}

	        if (bprintstatebyline)
	        {
	           //creates a file which has the state with the maximum posterior probability
	           String szmaxoutfilename = szoutputdir+"/STATEBYLINE/"+szprefix+"_"+szactualchrom+ChromHMM.SZSTATEBYLINEEXTENSION;
	           System.out.println("Writing to file "+szmaxoutfilename);
	           pwmax = new PrintWriter(szmaxoutfilename);
	           pwmax.println(cellSeq[nordered_nseq]+"\t"+szactualchrom);
	           pwmax.println("MaxState "+chorder);
		}

	        //PrintWriter pwbed = null;
	        if (bprintsegment)
	        {
	           //creates a file which has the maximum segmentation
	           //we only have one file per cell type here
	           pwbed = (PrintWriter) hmcellToFourColPW.get(cellSeq[nordered_nseq]);

	           if (pwbed == null)
	           {
		      //haven't seen this cell type
	              String szsegmentoutfilename = szoutputdir+"/" + szprefix+SZSEGMENTEXTENSION;

	              pwbed = new PrintWriter(szsegmentoutfilename);
		      System.out.println("Writing to file "+szsegmentoutfilename);
	              hmcellToFourColPW.put(cellSeq[nordered_nseq],pwbed);    
		   }
		}
	     }
	  }


          if (bscaleemissions)
	  {
	     for (int ni = 0; ni < emissionproducts.length; ni++)
	     {
	        //going through each combination of marks
		if (traindataObservedSeqFlags_nseq[ni])
		{
		   //this signature of marks is observed on the current chromosome so
		   //updating its emission probabilities
		   double[] emissionproducts_ni = emissionproducts[ni];
		   boolean[] traindataObservedValues_ni = traindataObservedValues[ni];
		   boolean[] traindataNotMissing_ni = traindataNotMissing[ni];

		   for (int ns = 0; ns < numstates; ns++)
		   {
		      emissionproducts_ni[ns] = 1;
		   }

		   for (int nmod = 0; nmod < numdatasets; nmod++)
		   {
		      for (int ns = 0; ns < numstates; ns++)
		      {
		         if (traindataNotMissing_ni[nmod])
			 {
			    //we are include this marks emission probability
			    if (traindataObservedValues_ni[nmod])
			    {
			       emissionproducts_ni[ns] *= emissionprobs[ns][nmod][1];
			    }
			    else
			    {
				emissionproducts_ni[ns] *= emissionprobs[ns][nmod][0];
			    }
			 }
			 // otherwise treated as missing omitting from product
		      }

		      double dmaxval = 0;
		      for (int ns = 0; ns < numstates; ns++)
		      {
		         if (emissionproducts_ni[ns] > dmaxval)
			 {
			    dmaxval = emissionproducts_ni[ns];
			 }
		      }

		      if (dmaxval <=0)
		      {
		         for (int ns = 0; ns < numstates; ns++)
	                 {
			    emissionproducts_ni[ns] = 1;
		         }
		      }
		      else
		      {
		         for (int ns = 0; ns < numstates; ns++)
	                 {
			    emissionproducts_ni[ns]/= dmaxval;
		         }
		      }
		   }
		}
	     }
	  }
	  else
	  {
	     for (int ni = 0; ni < emissionproducts.length; ni++)
             {
	        //going through each combination of marks
	        if (traindataObservedSeqFlags_nseq[ni])
	        {
	           //this signature of marks is observed on the current chromosome so
	           //updating its emission probabilities
	           double[] emissionproducts_ni = emissionproducts[ni];
	           boolean[] traindataObservedValues_ni = traindataObservedValues[ni];
	           boolean[] traindataNotMissing_ni = traindataNotMissing[ni];		  

		   boolean ballzero = true;

	          for (int ns = 0; ns < numstates; ns++)
	          {
	             double dproduct = 1;
	             double[][] emissionprobs_ni = emissionprobs[ns];

		     //going through all marks
		     for (int nmod = 0; nmod < numdatasets; nmod++)
	             {
		        if (traindataNotMissing_ni[nmod])
		        {
			   //we have observed the mark
		           if (traindataObservedValues_ni[nmod])
		           {
		              dproduct *= emissionprobs_ni[nmod][1];
			   }
		           else 
	                   {
		              dproduct *= emissionprobs_ni[nmod][0];
			   }
			}
		        // otherwise treated as missing omitting from product
		     }
	             emissionproducts_ni[ns] = dproduct;

		     if (dproduct >= EPSILONEMISSIONS)
		     {
		        ballzero = false;
		     }
		  }

		  if (ballzero)
	          {
	      	     for (int ns = 0; ns < numstates; ns++)
	             {
	       	        emissionproducts_ni[ns] = EPSILONEMISSIONS;
		     }
		  }
		}
	     }
	  }

	  //initial probability in state s is initial probability times emission probability at first position
          double[] alpha_nt = alpha[0];
	  double dscale = 0;
	  double[] emissionproducts_nobserveindex =emissionproducts[traindataObservedIndex_nseq[0]];
 	  for (int ns = 0; ns < numstates; ns++)
          {
	      alpha_nt[ns] = probinit[ns] * emissionproducts_nobserveindex[ns];
	      dscale += alpha_nt[ns];
	  }
	  scale[0] = dscale;

	  //alpha_t(s)=P(o_0,...,o_t,x_t=s|lambda)
          //converts the alpha terms to probabilities
          if (bscalebeta)
	  {
	     if (dscale == 0)
	     {
	        for (int ns = 0; ns < numstates; ns++)
	        {
	           //added
		   if ((alpha_nt[ns] < EPSILONSTATE)&&(emissionproducts_nobserveindex[ns]>0))//(ns < numstates-1))) 
	           {
	              alpha_nt[ns] = EPSILONSTATE;
		      dscale += EPSILONSTATE;
		   }
		}
	        scale[0] = dscale;
	     }

             for (int ns = 0; ns < numstates; ns++)
	     {
                alpha_nt[ns] /= dscale;

	        if ((alpha_nt[ns] < EPSILONSTATE)&&(emissionproducts_nobserveindex[ns]>0))//(ns < numstates-1))) 
       	        {
	           alpha_nt[ns] = EPSILONSTATE;
	        }
	     }
	  }
      	  else
	  {
             for (int ns = 0; ns < numstates; ns++)
	     {
                alpha_nt[ns] /= dscale;
	     }
	  }
	     //for (int ni = 0; ni < numstates; ni++)
	     //{
             //alpha_nt[ni] /= dscale;
	     //}
	
          //stores in coltransitionprobs the transpose of transitionprobs
          for (int ni = 0; ni < numstates; ni++)
          {
             double[] coltransitionprobs_ni = coltransitionprobs[ni];
             for (int nj = 0; nj < numstates; nj++)
	     {
	        coltransitionprobs_ni[nj] = transitionprobs[nj][ni];
	     }
	  }

          //forward step
          int numtime_nseq = numtime[nordered_nseq];

	  nprevlinecount = numtime_nseq;

          for (int nt = 1; nt < numtime_nseq; nt++)
          {
             //the actual observed combination at position t	        
	     double[] alpha_ntm1 = alpha[nt-1];
	     alpha_nt = alpha[nt];
	      
	     dscale = 0;
	     emissionproducts_nobserveindex = emissionproducts[traindataObservedIndex_nseq[nt]];
	     for (int ns = 0; ns < numstates; ns++)
	     {
	        //going through each state		   

	        int transitionprobsnumCol_ns = transitionprobsnumCol[ns];
	        int[] transitionprobsindexCol_ns = transitionprobsindexCol[ns];
	        double[] coltransitionprobs_ns = coltransitionprobs[ns];

	        double dtempsum = 0;
                if (transitionprobsnumCol_ns < nsparsecutoff)
	        {
		    //if it is sparse enough then it is worth the extra array indirection here
	           for (int nj = 0; nj < transitionprobsnumCol_ns; nj++)
	           {
	               //for each next state computing inner sum of all previous alpha and the transition probability
	               //for all non-zero transitions into the state
			int nmappedindex = transitionprobsindexCol_ns[nj];
			dtempsum += coltransitionprobs_ns[nmappedindex]*alpha_ntm1[nmappedindex];
		   }
		}
	        else
	        {
                   for (int nj = 0; nj < numstates; nj++)
	           {
	              //for each next state computing inner sum of all previous alpha and the transition probability
	              //for all transitions into the state
		      dtempsum += coltransitionprobs_ns[nj]*alpha_ntm1[nj];
		   }
		}

                //multiply the transition sum by the emission probability
	        double dalphaval = dtempsum*emissionproducts_nobserveindex[ns];
                alpha_nt[ns] = dalphaval;
	        dscale += dalphaval;
	     }

	      //rescaling alpha
              scale[nt] = dscale;
              //scale_t(s)=P(o_0,...,o_t|lambda) summed over all states

	      if (bscalebeta)
	      {
	         if (dscale == 0)
	         {
		    for (int ns = 0; ns < numstates; ns++)
		    {
		       //added
		       if ((alpha_nt[ns] < EPSILONSTATE)&&(emissionproducts_nobserveindex[ns]>0))//(ns < numstates-1))) 
		       {
		          alpha_nt[ns] = EPSILONSTATE;
			  dscale += EPSILONSTATE;
		       }
		    }
		    scale[nt] = dscale;
		 }

	         for (int ns = 0; ns < numstates; ns++)
	         {
                    alpha_nt[ns] /= dscale;

		    if ((alpha_nt[ns] < EPSILONSTATE)&&(emissionproducts_nobserveindex[ns]>0))//(ns < numstates-1))) 
		    {
		       alpha_nt[ns] = EPSILONSTATE;
		    }
		 }
	      }
      	      else
	      {
	         for (int ns = 0; ns < numstates; ns++)
	         {
                    alpha_nt[ns] /= dscale;
		 }
	      }
	     // for (int ns = 0; ns < numstates; ns++)
             // {
	     //	  alpha_nt[ns] /= dscale;
	     //}      	       
	  }
	    
          //backward step
          //beta_t(s)=P(o_t+1,...,o_T|x_t=s,lambda)
          int nlastindex = numtime_nseq-1;
	  double dinitval;
	  if (bscalebeta)
	  {
	     dinitval = 1.0/numstates;
	  }
	  else
	  {
             dinitval = 1.0/scale[nlastindex];
	  }

          for (int ns = 0; ns < numstates; ns++)
	  {
             beta_ntp1[ns] = dinitval;
	  }
          //double dinitval = 1.0/scale[nlastindex];
          //for (int ns = 0; ns < numstates; ns++)
          //{
          //    beta_ntp1[ns] = dinitval;
	  //}
	
	  int nmappedindexouter;
 
	  double ddenom = 0;	      

          //gamma_nt - P(x=S| o_0,...,o_t)
          //P(o_t+1,...,o_T|x_t=s,lambda) * P(o_0,...,o_t,x_t=s|lambda)
	  double[] gamma_nt = gamma[nlastindex]; 
          for (int ns = 0; ns < gamma_nt.length; ns++)
          {
	      double dval = alpha[nlastindex][ns]*beta_ntp1[ns];
	      ddenom += dval;
	      gamma_nt[ns] = dval;
	  }

	  if (ddenom > 0)
	  {
             for (int ns = 0; ns < gamma_nt.length; ns++)
             {
	        gamma_nt[ns] /= ddenom;
	     }
	  }


          for (int nt = nlastindex - 1; nt >= 0; nt--)
          {
	      gamma_nt = gamma[nt];
	      int ntp1 = (nt+1);
		   
	      double[] emissionproducts_ncombo_ntp1 = emissionproducts[traindataObservedIndex_nseq[ntp1]];		
	      double dsumbeta = 0;
	      double dscale_nt = scale[nt];

	      for (int ns = 0; ns < numstates; ns++)
              {
		  tempproductbetaemiss[ns] = beta_ntp1[ns]*emissionproducts_ncombo_ntp1[ns];
	      }

	      if (bscaleemissions)
	      {
		  //adding here to help numerical stability
		  double dmaxval = 0;
		  for (int ns = 0; ns < numstates; ns++)
		  {
		      if (tempproductbetaemiss[ns] > dmaxval)
		      {
			  dmaxval = tempproductbetaemiss[ns];
		      }
		  }
                  for (int ns = 0; ns < numstates; ns++)
		  {
		      tempproductbetaemiss[ns] = tempproductbetaemiss[ns]/dmaxval;
		  }
	      }

	      //double dscaleinv = 1.0/scale[nt];
              //scale_t(s)=P(o_0,...,o_t|lambda) summed over all states
	      for (int ni = 0; ni < numstates; ni++)
	      {
		  double dtempsum = 0;
		  int[] transitionprobsindex_ni =  transitionprobsindex[ni];
		  double[] transitionprobs_ni = transitionprobs[ni];
		  int transitionprobsnum_ni = transitionprobsnum[ni];

                  if (transitionprobsnum_ni < nsparsecutoff)
	          {
		    //if it is sparse enough then it is worth the extra array indirection here
	             for (int nj = 0; nj < transitionprobsnum_ni; nj++)
	             {
	                //for each state summing over transition probability to state j, emission probablity in j at next step
	                //and probability of observing the remaining sequence
		        nmappedindexouter = transitionprobsindex_ni[nj];
		        dtempsum += transitionprobs_ni[nmappedindexouter]*tempproductbetaemiss[nmappedindexouter];			
		     }
		  }
	          else
	          {
                     for (int nj = 0; nj < numstates; nj++)
	             {
	                //for each state summing over transition probability to state j, emission probablity in j at next step
	                //and probability of observing the remaining sequence
		        dtempsum += transitionprobs_ni[nj]*tempproductbetaemiss[nj];
		     }
		  }

		  if (bscalebeta)
		  {
		     beta_nt[ni] = dtempsum;
		     dsumbeta += dtempsum;
		  }
		  else
		  {
		     double dratio = dtempsum/dscale_nt;
		     if (dratio > Double.MAX_VALUE)
		     {
		        beta_nt[ni] = Double.MAX_VALUE;//dtempsum/dscale_nt;
		     }
		     else
		     {
		         beta_nt[ni] = dratio;
		     }
		  }
		  //double dratio = dtempsum/dscale_nt;
		  //if (dratio > Double.MAX_VALUE)
		  //{
		  //    beta_nt[ni] = Double.MAX_VALUE;
		  //}
		  //else
		  //{
		  //    beta_nt[ni] = dratio;//dtempsum/dscale_nt;
		  //}
	      }

	      if (bscalebeta)
	      {
                 for (int ni = 0; ni < numstates; ni++)
	         {
	            beta_nt[ni]/= dsumbeta;

		    if (beta_nt[ni] < EPSILONSTATE)//&&(!bdummy))// || (ni < numstates-1))) 
		    {
		       beta_nt[ni] = EPSILONSTATE;
		    }
		 }
	      }		

	      ddenom = 0;		
	      alpha_nt = alpha[nt];

	      //gamma_nt - P(x=S| o_0,...,o_t)
              //P(o_t+1,...,o_T|x_t=s,lambda) * P(o_0,...,o_t,xt=s|lambda)

	      for (int ns = 0; ns < gamma_nt.length; ns++)
              {
	         double dval = alpha_nt[ns]*beta_nt[ns];

		 ddenom += dval;
	         gamma_nt[ns] = dval;
	      }

	      if (ddenom > 0)
	      {
	         for (int ns = 0; ns < gamma_nt.length; ns++)
                 {
	            gamma_nt[ns]/=ddenom;       		   
		 }
	      }
	      beta_ntp1 = beta_nt;		
	  }


          gamma_nt = gamma[0];

	  double dmaxval = 0;
          int nmaxstate = 0;

	  int nt = 0;

	  if (bgzip)
	  {
	     //handling the first line
	     if (bnewfile)
	     {
                for (int ns = 0; ns < gamma_nt.length; ns++)
                {	   	
	           int nmappedstate = stateordering[ns]; //maps new state to old
	           double dprob = gamma_nt[nmappedstate];
	           if (bprintposterior)
	           {
		      String szout;
	              if (ns > 0)
	              {
		         szout = "\t"+nf.format(dprob);
		         //print with tab if not the first
                         //pwprobs.print("\t"+nf.format(dprob));
		      }
                      else
                      {
		         szout = nf.format(dprob);
                         //pwprobs.print(nf.format(dprob));
		      }

		      byte[] btformat = szout.getBytes();
		      pwprobszip.write(btformat,0,btformat.length);
		   }

	           if (dprob > dmaxval)
	           {
		      //best one found so far 
	              dmaxval = dprob;
	              nmaxstate = ns;
		   }
		}

	        if (bprintposterior)
	        {
		   String szout = "\n";
	  	   byte[] btformat = szout.getBytes();
		   pwprobszip.write(btformat,0,btformat.length);
	           //pwprobs.println();
		}

	        if (bprintstatebyline)
	        {
	           String szout = ""+(nmaxstate+1)+"\n";
	           byte[] btformat = szout.getBytes();
		   pwmaxzip.write(btformat,0,btformat.length);
	           //pwmax.println(""+(nmaxstate+1));
		}

	        //this contains the best state of the previous interval
	        nmaxstateprev = nmaxstate;  	
		nt = 1;
		nstart = 0;
	     }


             for ( ; nt < numtime_nseq; nt++)
             { 
                gamma_nt = gamma[nt];

	        dmaxval = 0;
	        nmaxstate = 0;
                for (int ns = 0; ns < gamma_nt.length; ns++)
                {	   	
		   int nmappedstate = stateordering[ns]; //maps new state to old
		   double dprob = gamma_nt[nmappedstate];
 	           if (bprintposterior)
	           {
		      String szout;
	              if (ns > 0)
	              {
			 szout = "\t"+nf.format(dprob);
		         //print with tab the first time
                         //pwprobs.print("\t"+nf.format(dprob));
		      }
	              else
	              {
			 szout = nf.format(dprob);
                         //pwprobs.print(nf.format(dprob));
		      }
		      byte[] btformat = szout.getBytes();
		      pwprobszip.write(btformat,0,btformat.length);
		   }

	           if (dprob > dmaxval)
	           {
	              dmaxval = dprob;
	              nmaxstate = ns;
		   }
		}

	        if (bprintposterior)
	        {
		    //pwprobs.println();
		   String szout = "\n";
		   byte[] btformat = szout.getBytes();
		   pwprobszip.write(btformat,0,btformat.length);
		}

                if (bprintstatebyline)
	        {
		   String szout = (nmaxstate+1)+"\n";
		   byte[] btformat = szout.getBytes();
		   pwmaxzip.write(btformat,0,btformat.length);
		   //pwmax.println(""+(nmaxstate+1));  	
		}	     
	     
	        if (bprintsegment&&(nmaxstateprev != nmaxstate))
	        {
		   //print out last segment we are done with
		    //String szout = chromSeq[nordered_nseq]+"\t"+((nstart+noffset)*nbinsize)+"\t"+((nt+noffset)*nbinsize)+"\t"+chorder+(nmaxstateprev+1) +"\n";
		   String szout = szactualchrom+"\t"+(nstart*nbinsize)+"\t"+((nt+noffset)*nbinsize)+"\t"+chorder+(nmaxstateprev+1) +"\n";
		    //pwbed.println(chromSeq[nordered_nseq]+"\t"+(nstart*nbinsize)+"\t"+(nt*nbinsize)+"\t"+chorder+(nmaxstateprev+1));	
                   byte[] btformat = szout.getBytes();
		   pwbedzip.write(btformat,0,btformat.length);
	 
		   //start a new segment now
		   nstart = nt+noffset;
		   nmaxstateprev = nmaxstate;
		}
	     }

	     if (bclosefile)
	     {
	        if (bprintsegment)
	        {
	           int nlastcoordinate;
	           Integer objMaxCoord = null;
	           if (hmMaxCoord != null)
	           {
		      //objMaxCoord = ((Integer) hmMaxCoord.get(chromSeq[nordered_nseq]));
		      objMaxCoord = ((Integer) hmMaxCoord.get(szactualchrom));
	           }

	           if (objMaxCoord != null)
	           {
		      nlastcoordinate = Math.min((numtime_nseq+noffset)*nbinsize,((Integer) objMaxCoord).intValue());

	  	      if (nstart*nbinsize >= nlastcoordinate)
		      {
			  //v1.26
			  throw new IllegalArgumentException("Start of interval "+(nstart*nbinsize+1)+" (1-based) is past the end of the chromosome "+nlastcoordinate+
                                                             "; Check if chromosome length file given to LearnModel is the same used for the binarization");
		      }
	           }
	           else
	           {
		      nlastcoordinate = (numtime_nseq+noffset)*nbinsize;
	           }
	           //pwbed.println(chromSeq[nordered_nseq]+"\t"+(nstart*nbinsize)+"\t"+nlastcoordinate+"\t"+chorder+(nmaxstateprev+1));
		   
		   String szout = szactualchrom+"\t"+(nstart*nbinsize)+"\t"+nlastcoordinate+"\t"+chorder+(nmaxstateprev+1)+"\n";
		   byte[] btformat = szout.getBytes();
		   pwbedzip.write(btformat,0,btformat.length);
	        }
	  
	        //close out the max state file if that was requested
	        if (bprintstatebyline)
	        {
	           pwmaxzip.finish();
	           pwmaxzip.close();
	        }
	        //close out the posterior state file if that was requested
	        if (bprintposterior)
	        {
	           pwprobszip.finish();
	           pwprobszip.close();    
	        }
	     }	        
	  }
	  else
	  {
	     if (bnewfile)
	     {
	        //handling the first line
                for (int ns = 0; ns < gamma_nt.length; ns++)
                {	   	
	           int nmappedstate = stateordering[ns]; //maps new state to old
	           double dprob = gamma_nt[nmappedstate];
	           if (bprintposterior)
	           {
	              if (ns > 0)
	              {
		         //print with tab if not the first
                         pwprobs.print("\t"+nf.format(dprob));
	              }
                      else
                      {
                         pwprobs.print(nf.format(dprob));
		      }
		   }

	           if (dprob > dmaxval)
	           {
		      //best one found so far 
	              dmaxval = dprob;
	              nmaxstate = ns;
		   }
	        }

	        if (bprintposterior)
	        {
	           pwprobs.println();
	        }

	        if (bprintstatebyline)
	        {
	           pwmax.println(""+(nmaxstate+1));
	        }

	        //this contains the best state of the previous interval
	        nmaxstateprev = nmaxstate;  	
		nt = 1;
		nstart = 0;
	     }

             for ( ; nt < numtime_nseq; nt++)
             { 
                gamma_nt = gamma[nt];

	        dmaxval = 0;
	        nmaxstate = 0;
                for (int ns = 0; ns < gamma_nt.length; ns++)
                {	   	
		   int nmappedstate = stateordering[ns]; //maps new state to old
		   double dprob = gamma_nt[nmappedstate];
 	           if (bprintposterior)
	           {
	              if (ns > 0)
	              {
		         //print with tab the first time
                         pwprobs.print("\t"+nf.format(dprob));
		      }
	              else
	              {
                         pwprobs.print(nf.format(dprob));
		      }
		   }

	           if (dprob > dmaxval)
	           {
	              dmaxval = dprob;
	              nmaxstate = ns;
		   }
		}

	        if (bprintposterior)
	        {
	           pwprobs.println();
		}

                if (bprintstatebyline)
	        {
		   pwmax.println(""+(nmaxstate+1));  	
		}	     
	     
	        if (bprintsegment&&(nmaxstateprev != nmaxstate))
	        {
		   //print out last segment we are done with
		    pwbed.println(szactualchrom+"\t"+(nstart*nbinsize)+"\t"+((nt+noffset)*nbinsize)+"\t"+chorder+(nmaxstateprev+1));		 
		   //start a new segment now
		   nstart = nt + noffset;
		   nmaxstateprev = nmaxstate;
		}
	     }

	     if (bclosefile)
	     {
	        if (bprintsegment)
	        {
	           int nlastcoordinate;
	           Integer objMaxCoord = null;
	           if (hmMaxCoord != null)
	           {
		      objMaxCoord = ((Integer) hmMaxCoord.get(szactualchrom));
	           }

	           if (objMaxCoord != null)
	           {
		      nlastcoordinate = Math.min((numtime_nseq+noffset)*nbinsize,((Integer) objMaxCoord).intValue());

	  	      if (nstart*nbinsize >= nlastcoordinate)
		      {
			  //v1.26
			  throw new IllegalArgumentException("Start of interval "+(nstart*nbinsize+1)+" (1-based) is past the end of the chromosome "+nlastcoordinate+
                                                             "; Check if chromosome length file given to LearnModel is the same used for the binarization");
		      }
	           }
	           else
	           {
		      nlastcoordinate = (numtime_nseq+noffset)*nbinsize;
	           }

	           pwbed.println(szactualchrom+"\t"+(nstart*nbinsize)+"\t"+nlastcoordinate+"\t"+chorder+(nmaxstateprev+1));
	        }
	  
	        //close out the max state file if that was requested
	        if (bprintstatebyline)
	        {
	           pwmax.close();
	        }
	        //close out the posterior state file if that was requested
	        if (bprintposterior)
	        {
	           pwprobs.close();    
		}
	     }	  
	  }
       }

       //if segment print was requested then we are going to go close those printwriters
       if (bprintsegment)
       {
          Iterator itr =  hmcellToFourColPW.values().iterator();
	  if (bgzip)
	  {
             while (itr.hasNext())
             {
	        GZIPOutputStream pwzip = (GZIPOutputStream) itr.next();
       	        pwzip.finish();
	        pwzip.close();
	     }
	  }
	  else 
          {
             //if segment print was requested then we are going to go close those printwriters
             while (itr.hasNext())
             {
	        PrintWriter pw = (PrintWriter) itr.next();
	        pw.close();
	     }	     
	  }	 
       }
    }


    ///////////////////////////////////////////////////////////////////////////////////////////

    /**
     * API call to compute max state at position npos (0-based) given input data
     */
    public int getMaxStateAtPos(int[][] data, int npos) throws IOException
    {

       if ((npos < 0) || (npos >= data.length))
       {
	   throw new IllegalArgumentException(npos+" is not in range of data"); 
       }

       //number of non-zero transition required to be less than this at the more stringent cutoff 
       //for trying to exploit sparsity in the transition matrix for efficiency gains
       int nsparsecutoff = (int) (numstates * ChromHMM.SPARSECUTOFFRATIO);

       //stores temporary product terms
       double[] tempproductbetaemiss = new double[numstates];

       //This stores the alpha values at each time point and number of states
       double[][] alpha = new double[data.length][numstates];

       //Temporary storage of the gamma's for each state
       //double[][] gamma = new double[data.length][numstates];

       //Temporary storage of the beta values for each state
       double[] beta_nt = new double[numstates];

       //Temporary storage of the beta values for each state at the next time point
       double[] beta_ntp1 = new double[numstates];

       //stores the scaling value for each time point
       double[] scale = new double[data.length];

       //stores the transition probabilities for each column
       double[][] coltransitionprobs = new double[numstates][numstates];
       

       //numdatasets is the number of marks we are integrating
       int numdatasets = data[0].length; 

       String[] szdataA = new String[data.length];
       for (int nrow = 0; nrow < data.length; nrow++)
       {
          StringBuffer sb = new StringBuffer();
		
	  int[] data_nrow = data[nrow];
	  for (int ncol = 0; ncol < numdatasets; ncol++)
	  {
	      int nval = data_nrow[ncol];
		    
	      if (nval == 0)
	      {
	         sb.append("0");
	      }
	      else if (nval == 1)
	      {
	         sb.append("1");
	      }
	      else if (nval == 2)
	      {
	         //this means missing
		 sb.append("2");
	      }
	      else
	      {
	         throw new IllegalArgumentException("Unrecognized value "+nval+" in input data");
	      }
	  }
	  szdataA[nrow] = sb.toString();
       }

       traindataObservedIndex = new int[1][data.length];
       int[] traindataObservedIndex_nseq = traindataObservedIndex[0];

       HashMap hmObserved = new HashMap();
       int nobserved = 0;

       for (int nrow = 0; nrow < szdataA.length; nrow++)
       {
          BigInteger theBigInteger = new BigInteger(szdataA[nrow],3);
	  ObservedRec theObservedRec  = (ObservedRec) hmObserved.get(theBigInteger);

	  if (theObservedRec == null)
	  {
	     //System.out.println(szmappingbyte.length());
	     //storing a mapping from observed byte string to an integer index in alFlags and alObserved
	     hmObserved.put(theBigInteger, new ObservedRec(nobserved,null));

	     //saving this observed index
	     traindataObservedIndex_nseq[nrow] = nobserved;

	     //increments the number of observed combinations of marks
	     nobserved++;
	  }
	  else
	  {
	     //storing the index of the flags associated with this row 
	     traindataObservedIndex_nseq[nrow] = theObservedRec.nobserved;		
	  }
       }
	    
       //saving the mapping of signatures and chromsome observed on

       //stores whether there is a present call at each location
       traindataObservedValues = new boolean[nobserved][numdatasets];
       //stores the emission probability for the i^th combination of marks in the j^th state
       double[][] emissionproducts = new double[traindataObservedValues.length][numstates];


       //stores whether the mark is not considered missing
       traindataNotMissing = new boolean[nobserved][numdatasets];

       Iterator hmObservedIterator = hmObserved.entrySet().iterator();
       while (hmObservedIterator.hasNext())
       {
          Map.Entry pairs = (Map.Entry) hmObservedIterator.next();
	  BigInteger theBigInteger = (BigInteger) pairs.getKey();
	  String szmapping = theBigInteger.toString(3);  //getting back the mapping string

	  ObservedRec theObservedRec = (ObservedRec) pairs.getValue();
	  int ncurrindex = theObservedRec.nobserved;//this is an index on which obervation combination it is

	  boolean[] traindataObservedValues_ncurrindex = traindataObservedValues[ncurrindex];
	  boolean[] traindataNotMissing_ncurrindex = traindataNotMissing[ncurrindex]; 
	    
	  //if the mapping string is less than the number of data sets then 
	  //there are leading 0's will set for leading 0's not missing and absent
	  int numch = szmapping.length();
	  int numleading0 = numdatasets - numch;
	  for (int nj = 0; nj < numleading0; nj++)
	  {
             traindataObservedValues_ncurrindex[nj] = false;
	     traindataNotMissing_ncurrindex[nj] = true;
	  }

	  int nmappedindex = numleading0; //starting from the leading 0 position
	  for (int nj = 0; nj < numch; nj++)
	  {
	     char ch = szmapping.charAt(nj);

	     if (ch == '0')
	     {
	        traindataObservedValues_ncurrindex[nmappedindex] = false;
	        traindataNotMissing_ncurrindex[nmappedindex] = true;
	     }
	     else if (ch=='1')
	     {
	        traindataObservedValues_ncurrindex[nmappedindex] = true;
		traindataNotMissing_ncurrindex[nmappedindex] = true;
	     }
	     else
	     {
	        //missing data
		traindataObservedValues_ncurrindex[nmappedindex] = false;
	        traindataNotMissing_ncurrindex[nmappedindex] = false;
	     }
	     nmappedindex++;
	  }
       }



       if (bscaleemissions)
       {
          for (int ni = 0; ni < emissionproducts.length; ni++)
	  {
	     //this signature of marks is observed on the current chromosome so
	     //updating its emission probabilities
	     double[] emissionproducts_ni = emissionproducts[ni];
	     boolean[] traindataObservedValues_ni = traindataObservedValues[ni];
	     boolean[] traindataNotMissing_ni = traindataNotMissing[ni];

	     for (int ns = 0; ns < numstates; ns++)
	     {
	        emissionproducts_ni[ns] = 1;
	     }

	     for (int nmod = 0; nmod < numdatasets; nmod++)
	     {
	        for (int ns = 0; ns < numstates; ns++)
	        {
	           if (traindataNotMissing_ni[nmod])
	           {
		      //we are include this marks emission probability
		      if (traindataObservedValues_ni[nmod])
		      {
		         emissionproducts_ni[ns] *= emissionprobs[ns][nmod][1];
		      }
		      else
		      {
		         emissionproducts_ni[ns] *= emissionprobs[ns][nmod][0];
		      }
		   }
			 // otherwise treated as missing omitting from product
		}

		double dmaxval = 0;
		for (int ns = 0; ns < numstates; ns++)
	        {
	           if (emissionproducts_ni[ns] > dmaxval)
		   {
	              dmaxval = emissionproducts_ni[ns];
		   }
		}

		if (dmaxval <=0)
	        {
	           for (int ns = 0; ns < numstates; ns++)
	           {
	              emissionproducts_ni[ns] = 1;
		   }
		}
	        else
	        {
	           for (int ns = 0; ns < numstates; ns++)
	           {
	              emissionproducts_ni[ns]/= dmaxval;
		   }
		}
	     }
	  }
       }	  
       else
       {
          for (int ni = 0; ni < emissionproducts.length; ni++)
          {
	     //this signature of marks is observed on the current chromosome so
	     //updating its emission probabilities
	     double[] emissionproducts_ni = emissionproducts[ni];
	     boolean[] traindataObservedValues_ni = traindataObservedValues[ni];
	     boolean[] traindataNotMissing_ni = traindataNotMissing[ni];		  

	     boolean ballzero = true;

             for (int ns = 0; ns < numstates; ns++)
	     {
	        double dproduct = 1;
	        double[][] emissionprobs_ni = emissionprobs[ns];

		//going through all marks
	        for (int nmod = 0; nmod < numdatasets; nmod++)
	        {
	           if (traindataNotMissing_ni[nmod])
		   {
	              //we have observed the mark
		      if (traindataObservedValues_ni[nmod])
		      {
		         dproduct *= emissionprobs_ni[nmod][1];
		      }
		      else 
	              {
		         dproduct *= emissionprobs_ni[nmod][0];
		      }
		   }
		   // otherwise treated as missing omitting from product
		}
	        emissionproducts_ni[ns] = dproduct;

		if (dproduct >= EPSILONEMISSIONS)
	        {
	           ballzero = false;
	        }
	     }

	     if (ballzero)
	     {
	        for (int ns = 0; ns < numstates; ns++)
	        {
		   emissionproducts_ni[ns] = EPSILONEMISSIONS;
	        }
	     }
	  }
       }

       //initial probability in state s is initial probability times emission probability at first position
       double[] alpha_nt = alpha[0];
       double dscale = 0;
       double[] emissionproducts_nobserveindex =emissionproducts[traindataObservedIndex_nseq[0]];
       for (int ns = 0; ns < numstates; ns++)
       {
	   alpha_nt[ns] = probinit[ns] * emissionproducts_nobserveindex[ns];
	   dscale += alpha_nt[ns];
       }
       scale[0] = dscale;

       //alpha_t(s)=P(o_0,...,o_t,x_t=s|lambda)
       //converts the alpha terms to probabilities
       if (bscalebeta)
       {
          if (dscale == 0)
          {
	     for (int ns = 0; ns < numstates; ns++)
	     {
	        //added
	        if ((alpha_nt[ns] < EPSILONSTATE)&&(emissionproducts_nobserveindex[ns]>0))//(ns < numstates-1))) 
	        {
	           alpha_nt[ns] = EPSILONSTATE;
	       	   dscale += EPSILONSTATE;
		}
	     }
	     scale[0] = dscale;
	  }

          for (int ns = 0; ns < numstates; ns++)
	  {
	     alpha_nt[ns] /= dscale;

	      if ((alpha_nt[ns] < EPSILONSTATE)&&(emissionproducts_nobserveindex[ns]>0))//(ns < numstates-1))) 
       	      {
	         alpha_nt[ns] = EPSILONSTATE;
	      }
	  }
       }
       else
       {
          for (int ns = 0; ns < numstates; ns++)
	  {
             alpha_nt[ns] /= dscale;
          }
       }
	
       //stores in coltransitionprobs the transpose of transitionprobs
       for (int ni = 0; ni < numstates; ni++)
       {
          double[] coltransitionprobs_ni = coltransitionprobs[ni];
          for (int nj = 0; nj < numstates; nj++)
	  {
             coltransitionprobs_ni[nj] = transitionprobs[nj][ni];
	  }
       }

       //forward step
       int numtime_nseq = data.length;


       for (int nt = 1; nt < numtime_nseq; nt++)
       {
          //the actual observed combination at position t	        
	  double[] alpha_ntm1 = alpha[nt-1];
	  alpha_nt = alpha[nt];
	      
	  dscale = 0;
	  emissionproducts_nobserveindex = emissionproducts[traindataObservedIndex_nseq[nt]];
	  for (int ns = 0; ns < numstates; ns++)
	  {
	     //going through each state		   
	     int transitionprobsnumCol_ns = transitionprobsnumCol[ns];
	     int[] transitionprobsindexCol_ns = transitionprobsindexCol[ns];
	     double[] coltransitionprobs_ns = coltransitionprobs[ns];

	     double dtempsum = 0;
             if (transitionprobsnumCol_ns < nsparsecutoff)
	     {
	        //if it is sparse enough then it is worth the extra array indirection here
	        for (int nj = 0; nj < transitionprobsnumCol_ns; nj++)
	        {
	           //for each next state computing inner sum of all previous alpha and the transition probability
	           //for all non-zero transitions into the state
		   int nmappedindex = transitionprobsindexCol_ns[nj];
		   dtempsum += coltransitionprobs_ns[nmappedindex]*alpha_ntm1[nmappedindex];
		}
	     }
	     else
	     {
                for (int nj = 0; nj < numstates; nj++)
	        {
	           //for each next state computing inner sum of all previous alpha and the transition probability
	           //for all transitions into the state
	           dtempsum += coltransitionprobs_ns[nj]*alpha_ntm1[nj];
		}
	     }

             //multiply the transition sum by the emission probability
	     double dalphaval = dtempsum*emissionproducts_nobserveindex[ns];
             alpha_nt[ns] = dalphaval;
	     dscale += dalphaval;
	  }

	  //rescaling alpha
          scale[nt] = dscale;
          //scale_t(s)=P(o_0,...,o_t|lambda) summed over all states

	  if (bscalebeta)
          {
	     if (dscale == 0)
	     {
	        for (int ns = 0; ns < numstates; ns++)
	        {
	           //added
	           if ((alpha_nt[ns] < EPSILONSTATE)&&(emissionproducts_nobserveindex[ns]>0))//(ns < numstates-1))) 
	           {
	              alpha_nt[ns] = EPSILONSTATE;
		      dscale += EPSILONSTATE;
		   }
		}
	        scale[nt] = dscale;
	     }

             for (int ns = 0; ns < numstates; ns++)
	     {
                alpha_nt[ns] /= dscale;

		if ((alpha_nt[ns] < EPSILONSTATE)&&(emissionproducts_nobserveindex[ns]>0))//(ns < numstates-1))) 
	        {
	           alpha_nt[ns] = EPSILONSTATE;
		}
	     }
	  }
      	  else
	  {
             for (int ns = 0; ns < numstates; ns++)
	     {
                alpha_nt[ns] /= dscale;
	     }
	  }
	     // for (int ns = 0; ns < numstates; ns++)
             // {
	     //	  alpha_nt[ns] /= dscale;
	     //}      	       
       }
	    
       //backward step
       //beta_t(s)=P(o_t+1,...,o_T|x_t=s,lambda)
       int nlastindex = numtime_nseq-1;
       double dinitval;
       if (bscalebeta)
       {
          dinitval = 1.0/numstates;
       }
       else
       {
         dinitval = 1.0/scale[nlastindex];
       }

       for (int ns = 0; ns < numstates; ns++)
       {
          beta_ntp1[ns] = dinitval;
       }

	
       int nmappedindexouter;


       //gamma_nt - P(x=S| o_0,...,o_t)
       //P(o_t+1,...,o_T|x_t=s,lambda) * P(o_0,...,o_t,x_t=s|lambda)
       //double[] gamma_nt = gamma[nlastindex]; 
       //for (int ns = 0; ns < gamma_nt.length; ns++)
       //{
       //   double dval = alpha[nlastindex][ns]*beta_ntp1[ns];
       //   ddenom += dval;
       //   gamma_nt[ns] = dval;
       //}

       //for (int ns = 0; ns < gamma_nt.length; ns++)
       //{
       //   gamma_nt[ns] /= ddenom;
       //}
 

       for (int nt = nlastindex - 1; nt >= npos; nt--) //only goes to npos not 0
       {
	   //gamma_nt = gamma[nt];
	  int ntp1 = (nt+1);
		   
	  double[] emissionproducts_ncombo_ntp1 = emissionproducts[traindataObservedIndex_nseq[ntp1]];		
          double dsumbeta = 0;
	  double dscale_nt = scale[nt];

	  for (int ns = 0; ns < numstates; ns++)
          {
             tempproductbetaemiss[ns] = beta_ntp1[ns]*emissionproducts_ncombo_ntp1[ns];
	  }

	  if (bscaleemissions)
          {
	     //adding here to help numerical stability
	     double dmaxval = 0;
	     for (int ns = 0; ns < numstates; ns++)
	     {
                if (tempproductbetaemiss[ns] > dmaxval)
	        {
	       	   dmaxval = tempproductbetaemiss[ns];
	        }
	     }
             for (int ns = 0; ns < numstates; ns++)
	     {
                tempproductbetaemiss[ns] = tempproductbetaemiss[ns]/dmaxval;
	     }
	  }

	  //double dscaleinv = 1.0/scale[nt];
          //scale_t(s)=P(o_0,...,o_t|lambda) summed over all states
	  for (int ni = 0; ni < numstates; ni++)
	  {
	     double dtempsum = 0;
	     int[] transitionprobsindex_ni =  transitionprobsindex[ni];
	     double[] transitionprobs_ni = transitionprobs[ni];
	     int transitionprobsnum_ni = transitionprobsnum[ni];

             if (transitionprobsnum_ni < nsparsecutoff)
	     {
	        //if it is sparse enough then it is worth the extra array indirection here
	        for (int nj = 0; nj < transitionprobsnum_ni; nj++)
	        {
	           //for each state summing over transition probability to state j, emission probablity in j at next step
	           //and probability of observing the remaining sequence
		   nmappedindexouter = transitionprobsindex_ni[nj];
	           dtempsum += transitionprobs_ni[nmappedindexouter]*tempproductbetaemiss[nmappedindexouter];			
		}
	     }
	     else
	     {
                for (int nj = 0; nj < numstates; nj++)
	        {
	           //for each state summing over transition probability to state j, emission probablity in j at next step
	           //and probability of observing the remaining sequence
		   dtempsum += transitionprobs_ni[nj]*tempproductbetaemiss[nj];
		}
	     }

	     if (bscalebeta)
	     {
	        beta_nt[ni] = dtempsum;
       	        dsumbeta += dtempsum;
	     }
	     else
	     {
	        double dratio = dtempsum/dscale_nt;
		if (dratio > Double.MAX_VALUE)
	        {
	           beta_nt[ni] = Double.MAX_VALUE;//dtempsum/dscale_nt;
	        }
	        else
	        {
	           beta_nt[ni] = dratio;
      	        }
	     }
		  //double dratio = dtempsum/dscale_nt;
		  //if (dratio > Double.MAX_VALUE)
		  //{
		  //    beta_nt[ni] = Double.MAX_VALUE;
		  //}
		  //else
		  //{
		  //    beta_nt[ni] = dratio;//dtempsum/dscale_nt;
		  //}
	  }

	  if (bscalebeta)
	  {
             for (int ni = 0; ni < numstates; ni++)
	     {
	        beta_nt[ni]/= dsumbeta;

		if (beta_nt[ni] < EPSILONSTATE)//&&(!bdummy))// || (ni < numstates-1))) 
	        {
	           beta_nt[ni] = EPSILONSTATE;
	        }
	     }
	  }	

	  //ddenom = 0;		
          //alpha_nt = alpha[nt];
	  beta_ntp1 = beta_nt;
       }

       //gamma_nt - P(x=S| o_0,...,o_t)
       //P(o_t+1,...,o_T|x_t=s,lambda) * P(o_0,...,o_t,xt=s|lambda)

       //for (int ns = 0; ns < gamma_nt.length; ns++)
       //{
       //   double dval = alpha_nt[ns]*beta_nt[ns];
       //    ddenom += dval;
       //   gamma_nt[ns] = dval;
       //}

       //for (int ns = 0; ns < gamma_nt.length; ns++)
       //{
       //  gamma_nt[ns]/=ddenom;       		  
       //}
       //}

       double dmaxval = 0;
       int nmaxstate = 0;
       double[] alpha_npos = alpha[npos];
       //double[] beta_npos = beta[npos];
       //double[] gamma_npos = gamma[npos];


       for (int ns = 0; ns < alpha_npos.length; ns++)
       {
	  double dval = alpha_npos[ns]*beta_ntp1[ns];// gamma_npos[ns];
	  if (dval > dmaxval)
          {
             //best one found so far 
	     dmaxval = dval;
	     nmaxstate = ns;
          }
       }

       return (nmaxstate+1);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    /**
     * This is the core procedure for learning the parameters of the model
     */
    public void trainParametersWithLoad() throws IOException
    {

        NumberFormat nf3 = NumberFormat.getInstance();
        nf3.setMaximumFractionDigits(3);
	nf3.setGroupingUsed(false);
	nf3.setMinimumFractionDigits(3);

        NumberFormat nf1 = NumberFormat.getInstance();
	nf1.setMaximumFractionDigits(1);
	nf1.setMinimumFractionDigits(1);
	nf1.setGroupingUsed(false);

       int niteration = 1;

       boolean bconverged = false;

       double dzerotransitioncutoff = Math.pow(10,-nzerotransitionpower);

       //number of non-zero transition for the 
       int nsparsecutoff = (int) (numstates * ChromHMM.SPARSECUTOFFRATIO);

       //number of non-zero transition that need to be less than this at the looser cut-off
       int nsparsecutofflooser = (int) (numstates * ChromHMM.SPARSECUTOFFLOOSERRATIO);

       double dprevloglike;

       //stores the maximum number of locations in any sequence and in each sequence
       //int[] numtime = new int[traindataObservedIndex.length];
       int nmaxtime = 0;
       for (int nseq = 0; nseq < numtime.length; nseq++)
       {
	   //numtime[nseq] = traindataObservedIndex[nseq].length;
          if (numtime[nseq] > nmaxtime)
	  {
      	     nmaxtime = numtime[nseq];
	  }
       }

       if (ChromHMM.BVERBOSE)
       {
          System.out.println("Maximum number of locations\t"+nmaxtime);
       }

       //stores the emission probability for the i^th combination of marks in the j^th state
       double[][] emissionproducts = new double[nmaxtime][numstates]; //changed to maxtime here

       //stores the scale for the emission probability for the i^th combination of marks
       double[] emissionproducts_scale = null;

       if (bscaleemissions)
       {
          emissionproducts_scale = new double[nmaxtime];
       }

       //stores temporary product terms
       double[] tempproductbetaemiss = new double[numstates];

       //This stores the alpha values at each time point and number of states
       double[][] alpha = new double[nmaxtime][numstates];

       //Temporary storage of the gamma's for each state
       double[] gamma_nt = new double[numstates];

       //Temporary storage of the beta values for each state
       double[] beta_nt = new double[numstates];

       //Temporary storage of the beta values for each state at the next time point
       double[] beta_ntp1 = new double[numstates];

       //stores the scaling value for each time point
       double[] scale = new double[nmaxtime];

       //stores the transition probabilities for each column
       double[][] coltransitionprobs = new double[numstates][numstates];

       //stores the sufficient statistic for the initital probability in each state for the last visit
       double[][] gammainitstore = new double[chromfiles.length][numstates];//changed from traindataObservedIndex.length

       //stores the sufficient statistics for computing the transition probabilities cumulated for each iteration
       double[][][] sxistore = new double[chromfiles.length][numstates][numstates];

       //stores the sufficient statistic for computing the emission probabilities
       double[][][][] gammaksumstore = 
                 new double[chromfiles.length][numstates][numdatasets][numbuckets];

       //temporary storage in computation of sxi
       double[][] sumforsxi = new double[numstates][numstates];

       //stores the sum of the gamma values associated with each combination in each state
       double[][] gammaObservedSum = new double[nmaxtime][numstates];//changed from traindataObservedValues.length

       //stores the indicies of the data
       int[] traindataObservedIndex = new int[nmaxtime];

      //saving the mapping of signatures and chromsome observed on
      //stores whether there is a present call at each location
      boolean[][] traindataObservedValues = new boolean[nmaxtime][numdatasets];

      //stores whether the mark is not considered missing
      boolean[][] traindataNotMissing = new boolean[nmaxtime][numdatasets]; //usually nobserved

       int nelim = 0;
       long ltimeitr= System.currentTimeMillis();
       dprevloglike = Double.NEGATIVE_INFINITY;

       do
       {
          dloglike = 0;	  

          //for (int nseq = 0; nseq < traindataObservedIndex.length; nseq++)
	  for (int nseq = 0; nseq < chromfiles.length; nseq++)
          {
	      //going through each sequence


	      ///////////////////////////////////////////////////////////
	      //loading the data just for the current chromosome 

	      HashMap hmObserved = new HashMap(); //maps an observation string to an index and set of flags
 
    	      int nobserved = 0;

	      if (ChromHMM.BVERBOSE)
	      {
                 System.out.println("reading\t"+szinputdir+" "+chromfiles[nseq]);
	      }
	      BufferedReader br = Util.getBufferedReader(szinputdir+"/"+chromfiles[nseq]);
              String szLine; 
              br.readLine(); //first line tells cell type and chromosome
	      br.readLine();//flush mark header
	      ArrayList aldata = new ArrayList();
	      while ((szLine = br.readLine())!=null)
	      {
	         StringTokenizer st = new StringTokenizer(szLine,"\t ");
		 StringBuffer sb = new StringBuffer();
		
		 for (int ncol = 0; ncol < numdatasets; ncol++)
		 {
		    if (!st.hasMoreTokens())
		    {
	      	       throw new IllegalArgumentException("Found line without "+numdatasets+" values in file "+chromfiles[nseq]);
		    }

		    String sztoken = st.nextToken();//.trim();
		    
		    if (sztoken.equals("0"))
		    {
		       sb.append("0");
		    }
		    else if (sztoken.equals("1"))
		    {
		       sb.append("1");
		    }
		    else if (sztoken.equals("2"))
		    {
		       //this means missing
	               sb.append("2");
		    }
		    else
	            {
	               throw new IllegalArgumentException("Unrecognized value "+sztoken+" found in "+szinputdir+"/"+chromfiles[nseq]);
	            }
		 }
		 aldata.add(sb.toString());
	      }
	      br.close();
	      
	      int nsize = aldata.size();
	      //traindataObservedIndex[nfile] = new int[nsize];
	      //int[] traindataObservedIndex_nfile = traindataObservedIndex[nfile];

	      for (int nrow = 0; nrow < nsize; nrow++)
	      {
	         BigInteger theBigInteger = new BigInteger((String) aldata.get(nrow),3);
		 Integer theObservedInt  = (Integer) hmObserved.get(theBigInteger);
	         //boolean[] flagA;

	         if (theObservedInt == null)
       	         {
		    //System.out.println(szmappingbyte.length());
		    //storing a mapping from observed byte string to an integer index in alFlags and alObserved
		     hmObserved.put(theBigInteger, Integer.valueOf(nobserved));

		    //saving this observed index
		    traindataObservedIndex[nrow] = nobserved;

		    //increments the number of observed combinations of marks
		    nobserved++;
		 }
		 else
		 {
		    //storing the index of the flags associated with this row 
		     traindataObservedIndex[nrow] = ((Integer) theObservedInt).intValue();
		 }
	      }	  	 
	    




	     Iterator hmObservedIterator = hmObserved.entrySet().iterator();
	     while (hmObservedIterator.hasNext())
	     {
	        Map.Entry pairs = (Map.Entry) hmObservedIterator.next();
	        BigInteger theBigInteger = (BigInteger) pairs.getKey();
	        String szmapping = theBigInteger.toString(3);  //getting back the mapping string

	        //ObservedRec theObservedRec = (ObservedRec) pairs.getValue();
	        int ncurrindex = ((Integer) pairs.getValue()).intValue();// theObservedRec.nobserved;//this is an index on which obervation combination it is

	        boolean[] traindataObservedValues_ncurrindex = traindataObservedValues[ncurrindex];
	        boolean[] traindataNotMissing_ncurrindex = traindataNotMissing[ncurrindex]; 
	    
	        //if the mapping string is less than the number of data sets then 
	        //there are leading 0's will set for leading 0's not missing and absent
	        int numch = szmapping.length();
	        int numleading0 = numdatasets - numch;
	        for (int nj = 0; nj < numleading0; nj++)
	        {
		   traindataObservedValues_ncurrindex[nj] = false;
		   traindataNotMissing_ncurrindex[nj] = true;
		}

	        int nmappedindex = numleading0; //starting from the leading 0 position
	        for (int nj = 0; nj < numch; nj++)
	        {
	           char ch = szmapping.charAt(nj);

	           if (ch == '0')
	           {
		      traindataObservedValues_ncurrindex[nmappedindex] = false;
		      traindataNotMissing_ncurrindex[nmappedindex] = true;
		   }
	           else if (ch=='1')
	           {
		      traindataObservedValues_ncurrindex[nmappedindex] = true;
		      traindataNotMissing_ncurrindex[nmappedindex] = true;
		   }
	           else
	           {
		      //missing data
		      traindataObservedValues_ncurrindex[nmappedindex] = false;
		      traindataNotMissing_ncurrindex[nmappedindex] = false;
		   }
	           nmappedindex++;
		}
	     }      




	      ///////////////////////////////////////////////////////////
	      //int[] traindataObservedIndex_nseq = traindataObservedIndex[nseq];
	      //boolean[] traindataObservedSeqFlags_nseq = traindataObservedSeqFlags[nseq];

             double[][][] gammaksum_nseq = gammaksumstore[nseq];

	     for (int ns = 0; ns < gammaksum_nseq.length; ns++)
	     {
		 //resetting the gamma sufficient statistics in the current sequence
		double[][] gammaksum_nseq_ns = gammaksum_nseq[ns];
	        for (int nmark = 0; nmark < gammaksum_nseq_ns.length; nmark++)
	        {
		   for (int nbucket = 0; nbucket < numbuckets; nbucket++)
		   {
	              gammaksum_nseq_ns[nmark][nbucket] = 0;
	           }
	        }
	     }

	     double[][] sxi_nseq = sxistore[nseq];
	     for (int ni = 0; ni < sxi_nseq.length; ni++)
	     {
		 //reseeting the sxi sufficient statistics in the current sequence
		double[] sxi_nseq_ni = sxi_nseq[ni];
	        for (int nj = 0; nj < sxi_nseq_ni.length; nj++)
	        {
	           sxi_nseq_ni[nj] = 0;
	        }
	     }

	     //gammaObservedSum stores the weight for each combination of marks in each state
	     //for (int ncombo = 0; ncombo < gammaObservedSum.length; ncombo++)
	     for (int ncombo = 0; ncombo < nobserved; ncombo++)
	     {
		 //resetting that to 0
		double[] gammaObservedSum_ncombo = gammaObservedSum[ncombo];
	        for (int ns = 0; ns < gammaObservedSum_ncombo.length; ns++)
	        {
	           gammaObservedSum_ncombo[ns] = 0;
	        }
	     }


             if (bscaleemissions)
	     {
	        //for (int ni = 0; ni < emissionproducts.length; ni++)
		for (int ni = 0; ni < nobserved; ni++)
	        {
	           //going through each combination of marks
		   //if (traindataObservedSeqFlags_nseq[ni])
		   //{
		   //this signature of marks is observed on the current chromosome so
		   //updating its emission probabilities
		   double[] emissionproducts_ni = emissionproducts[ni];
		   boolean[] traindataObservedValues_ni = traindataObservedValues[ni];
		   boolean[] traindataNotMissing_ni = traindataNotMissing[ni];

		   for (int ns = 0; ns < numstates; ns++)
		   {
		       emissionproducts_ni[ns] = 1;
		   }
		   emissionproducts_scale[ni] = 0;

		   for (int nmod = 0; nmod < numdatasets; nmod++)
		   {
		      for (int ns = 0; ns < numstates; ns++)
		      {
		         if (traindataNotMissing_ni[nmod])
			 {
			    //we are include this marks emission probability
			    if (traindataObservedValues_ni[nmod])
			    {
			       emissionproducts_ni[ns] *= emissionprobs[ns][nmod][1];
			    }
			 }
					     // otherwise treated as missing omitting from product
		      }

		      double dmaxval = 0;
		      for (int ns = 0; ns < numstates; ns++)
		      {
		         if (emissionproducts_ni[ns] > dmaxval)
			 {
			     dmaxval = emissionproducts_ni[ns];
			 }
		      }

		      if (dmaxval <= 0)//EPSILONEMISSIONS)
		      {
		         for (int ns = 0; ns < numstates; ns++)
		         {
			    emissionproducts_ni[ns] = 1;
			 }
			 //approximating log-likelihood
                         emissionproducts_scale[ni] += Math.log(EPSILONEMISSIONS);
		      }
		      else
		      {
		         for (int ns = 0; ns < numstates; ns++)
		         {
			    emissionproducts_ni[ns]/= dmaxval;
		         }
		         emissionproducts_scale[ni] += Math.log(dmaxval);
		      }
		   }
		}
	     }
	     else
	     {
	        //for (int ni = 0; ni < emissionproducts.length; ni++)
	        for (int ni = 0; ni < nobserved; ni++)
	        {
	           //going through each combination of marks
	           //if (traindataObservedSeqFlags_nseq[ni])
		   //{
	           //this signature of marks is observed on the current chromosome so
       	           //updating its emission probabilities
       	           double[] emissionproducts_ni = emissionproducts[ni];
      	           boolean[] traindataObservedValues_ni = traindataObservedValues[ni];
       	           boolean[] traindataNotMissing_ni = traindataNotMissing[ni];		  

	           boolean ballzero = true;

	           for (int ns = 0; ns < numstates; ns++)
	           {
	              double dproduct = 1;
	              double[][] emissionprobs_ns = emissionprobs[ns];

		      for (int nmod = 0; nmod < numdatasets; nmod++)
	              {
		         if (traindataNotMissing_ni[nmod])
	      	         {
		            //we are include this marks emission probability
		            if (traindataObservedValues_ni[nmod])
		            {
		      	       //System.out.println("positive\t"+ns+"\t"+nmod+"\t1\t"+emissionprobs_ns[nmod][1]);
		               dproduct *= emissionprobs_ns[nmod][1];
		            }
		            else 
		            {
		      	       ///System.out.println("negative\t"+ns+"\t"+nmod+"\t0\t"+emissionprobs_ns[nmod][0]);
		               dproduct *= emissionprobs_ns[nmod][0];
			    }
			 }
			   // otherwise treated as missing omitting from product
		      }
		      //System.out.println(ns+"\t"+dproduct);
		      emissionproducts_ni[ns] = dproduct;

		      if (dproduct >= EPSILONEMISSIONS)
		      {
	                 ballzero = false;
		      }
		   }

		   if (ballzero)
	           {
	              for (int ns = 0; ns < numstates; ns++)
		      {
		         emissionproducts_ni[ns] = EPSILONEMISSIONS;
		      }
		   }
		}
	     }
	  
	     //initial probability in state s is initial probability times emission probability at first position
	     double[] alpha_nt = alpha[0];
	     double dscale = 0;
	     double[] emissionproducts_nobserveindex =emissionproducts[traindataObservedIndex[0]];
 	     for (int ns = 0; ns < numstates; ns++)
             {
	        alpha_nt[ns] = probinit[ns] * emissionproducts_nobserveindex[ns];
	        dscale += alpha_nt[ns];
	     }
	     scale[0] = dscale;

	     //alpha_t(s)=P(o_0,...,o_t,x_t=s|lambda)
	     //converts the alpha terms to probabilities
	     if (bscalebeta)
	     {
	        if (dscale == 0)
	        {
		   for (int ns = 0; ns < numstates; ns++)
		   {
		      //added
		      if ((alpha_nt[ns] < EPSILONSTATE)&&(emissionproducts_nobserveindex[ns]>0))//(ns < numstates-1))) 
		      {
		         alpha_nt[ns] = EPSILONSTATE;
		         dscale += EPSILONSTATE;
		      }
		   }
		   scale[0] = dscale;
		}

	        for (int ns = 0; ns < numstates; ns++)
	        {
                   alpha_nt[ns] /= dscale;

		   if ((alpha_nt[ns] < EPSILONSTATE)&&(emissionproducts_nobserveindex[ns]>0))//(ns < numstates-1))) 
		   {
		      alpha_nt[ns] = EPSILONSTATE;
		   }
		}
	     }
      	     else
	     {
	        for (int ns = 0; ns < numstates; ns++)
	        {
                   alpha_nt[ns] /= dscale;
		}
	     }
	     //for (int ni = 0; ni < numstates; ni++)
	     //{
	     //	 alpha_nt[ni] /= dscale;
	     //}
	     dloglike += Math.log(dscale); 

             if (bscaleemissions)
	     {
	        dloglike += emissionproducts_scale[traindataObservedIndex[0]];
	     }

	     //stores in coltransitionprobs the transpose of transitionprobs
	     for (int ni = 0; ni < numstates; ni++)
	     {
	        double[] coltransitionprobs_ni = coltransitionprobs[ni];
	        for (int nj = 0; nj < numstates; nj++)
	        {
	           coltransitionprobs_ni[nj] = transitionprobs[nj][ni];
	        }
	     }

	     //forward step
	     int numtime_nseq = numtime[nseq];
	     for (int nt = 1; nt < numtime_nseq; nt++)
	     {
	        //the actual observed combination at position t	        
	        double[] alpha_ntm1 = alpha[nt-1];
	        alpha_nt = alpha[nt];
	      
	        dscale = 0;
		emissionproducts_nobserveindex = emissionproducts[traindataObservedIndex[nt]];
	        for (int ns = 0; ns < numstates; ns++)
	        {
		   //stores the emission product for each location on the chromosome		   

		   int transitionprobsnumCol_ns = transitionprobsnumCol[ns];
		   int[] transitionprobsindexCol_ns = transitionprobsindexCol[ns];
	           double[] coltransitionprobs_ns = coltransitionprobs[ns];

	           double dtempsum = 0;
                   if (transitionprobsnumCol_ns < nsparsecutoff)
		   {
		       //number of transitions is sparse enough worth going through the extra redirection
	              for (int nj = 0; nj < transitionprobsnumCol_ns; nj++)
		      {
		         //for each next state computing inner sum of all previous alpha and the transition probability
		         //for all non-zero transitions into the state
		         int nmappedindex = transitionprobsindexCol_ns[nj];
		         dtempsum += coltransitionprobs_ns[nmappedindex]*alpha_ntm1[nmappedindex];
		      }
		   }
		   else
	           {
		       //avoid the redirect and multiply by 0
                      for (int nj = 0; nj < numstates; nj++)
		      {
		         dtempsum += coltransitionprobs_ns[nj]*alpha_ntm1[nj];
		      }
		   }
	           //multiply the transition sum by the emission probability
		   double dalphaval = dtempsum*emissionproducts_nobserveindex[ns];
	           alpha_nt[ns] = dalphaval;
		   //System.out.println(ns+"\t"+alpha_nt[ns]+"\t"+dtempsum+"\t"+emissionproducts_nobserveindex[ns]);
	           dscale += dalphaval;
	        }


		//rescaling alpha
	        scale[nt] = dscale;
                //scale_t(s)=P(o_0,...,o_t|lambda) summed over all states

	        if (bscalebeta)
	        {
	           if (dscale == 0)
	           {
		      for (int ns = 0; ns < numstates; ns++)
		      {
		         //added
		         if ((alpha_nt[ns] < EPSILONSTATE)&&(emissionproducts_nobserveindex[ns]>0))//(ns < numstates-1))) 
		         {
		            alpha_nt[ns] = EPSILONSTATE;
			    dscale += EPSILONSTATE;
		         }
		      }
		      scale[nt] = dscale;
		   }

	           for (int ns = 0; ns < numstates; ns++)
	           {
                      alpha_nt[ns] /= dscale;

		      if ((alpha_nt[ns] < EPSILONSTATE)&&(emissionproducts_nobserveindex[ns]>0))//(ns < numstates-1))) 
		      {
		         alpha_nt[ns] = EPSILONSTATE;
		      }
		   }
	        }
      	        else
	        {
	           for (int ns = 0; ns < numstates; ns++)
	           {
                      alpha_nt[ns] /= dscale;
		   }
	        }
	        //for (int ns = 0; ns < numstates; ns++)
	        //{
		//  alpha_nt[ns] /= dscale;
		//}

      	        dloglike += Math.log(dscale);

                if (bscaleemissions)
		{
		   dloglike += emissionproducts_scale[traindataObservedIndex[nt]];
	        }
	     }
	    
	     //backward step
	     //beta_t(s)=P(o_t+1,...,o_T|x_t=s,lambda)
	     int nlastindex = numtime_nseq-1;
	     double dinitval;
	     if (bscalebeta)
	     {
		dinitval = 1.0/numstates;
	     }
	     else
	     {
                dinitval = 1.0/scale[nlastindex];
	     }

             for (int ns = 0; ns < numstates; ns++)
	     {
	        beta_ntp1[ns] = dinitval;
	     }
	     //double dinitval = 1.0/scale[nlastindex];
             //for (int ns = 0; ns < numstates; ns++)
	     //{
	     //  beta_ntp1[ns] = dinitval;
	     //}

	     double ddenom = 0;	       

	     //gamma_nt - P(x=S| o_0,...,o_t)
	     //P(o_t+1,...,o_T|x_t=s,lambda) * P(o_0,...,o_t,xt=s|lambda)
	     alpha_nt = alpha[nlastindex];
	     for (int ns = 0; ns < gamma_nt.length; ns++)
	     {
	        double dval = alpha_nt[ns]*beta_ntp1[ns];
	        ddenom += dval;
	        gamma_nt[ns] = dval;
	     }

	     if (ddenom > 0)
	     {
                for (int ns = 0; ns < gamma_nt.length; ns++)
	        {
		   gamma_nt[ns] /=ddenom;
		}
	     }

	     double[] gammaObservedSum_combo_nt = gammaObservedSum[traindataObservedIndex[nlastindex]];
		
	     for (int ns = 0; ns < numstates; ns++)
	     { 
                //first sum gamma over all common signatures		     
	        //updates probability of observing the signature when in the state
	        gammaObservedSum_combo_nt[ns] += gamma_nt[ns];
	     }

	     for (int nt = nlastindex - 1; nt >= 0; nt--)
	     {
	        int ntp1 = (nt+1);
		   
	        double[] emissionproducts_combo_ntp1 = emissionproducts[traindataObservedIndex[ntp1]];

		for (int ns = 0; ns < numstates; ns++)
	        {
		    tempproductbetaemiss[ns] = beta_ntp1[ns]*emissionproducts_combo_ntp1[ns];
		}

	        if (bscaleemissions)
	        {
		   //adding here to help numerical stability
		   double dmaxval = 0;
		   for (int ns = 0; ns < numstates; ns++)
		   {
		      if (tempproductbetaemiss[ns] > dmaxval)
		      {
			  dmaxval = tempproductbetaemiss[ns];
		      }
		   }
                   for (int ns = 0; ns < numstates; ns++)
		   {
		      tempproductbetaemiss[ns] = tempproductbetaemiss[ns]/dmaxval;
		   }
	        }

		//double dscaleinv = 1.0/scale[nt];
		double dsumbeta = 0;
		double dscale_nt = scale[nt];
                //scale_t(s)=P(o_0,...,o_t|lambda) summed over all states

	        for (int ni = 0; ni < numstates; ni++)
		{
	           double dtempsum = 0;
	           int[] transitionprobsindex_ni =  transitionprobsindex[ni];
	           double[] transitionprobs_ni = transitionprobs[ni];
		   int transitionprobsnum_ni = transitionprobsnum[ni];

                   if (transitionprobsnum_ni < nsparsecutoff)
		   {
		       //sparse enought to pay the indirection penalty
		      for (int nj = 0; nj < transitionprobsnum_ni; nj++)
		      {
		         //for each state summing over transition probability to state j, emission probablity in j at next step
		         //and probability of observing the remaining sequence
			 int nmappedindexouter = transitionprobsindex_ni[nj];
		         dtempsum += transitionprobs_ni[nmappedindexouter]*tempproductbetaemiss[nmappedindexouter];
		      }
		   }
		   else
		   {
		       //not trying to exploit sparsity here
                      for (int nj = 0; nj < numstates; nj++)
		      {
		         //for each state summing over transition probability to state j, emission probablity in j at next step
		         //and probability of observing the remaining sequence
			 dtempsum += transitionprobs_ni[nj]*tempproductbetaemiss[nj];
		      }
		   }

		   if (bscalebeta)
		   {
		       beta_nt[ni] = dtempsum;
		       dsumbeta += dtempsum;
		   }
		   else
		   {
		      double dratio = dtempsum/dscale_nt;
		      if (dratio > Double.MAX_VALUE)
		      {
		         beta_nt[ni] = Double.MAX_VALUE;//dtempsum/dscale_nt;
		      }
		      else
		      {
		         beta_nt[ni] = dratio;
		      }
		   }
		   //double dratio = dtempsum/dscale_nt;
		   //if (dratio > Double.MAX_VALUE)
		   //{
		   //   beta_nt[ni] = Double.MAX_VALUE;
		   //}
		   //else
		   //{
		   //   beta_nt[ni] = dratio;
		   //}
		}		
		
	        if (bscalebeta)
	        {
                   for (int ni = 0; ni < numstates; ni++)
	           {
	               beta_nt[ni]/= dsumbeta;

		       if (beta_nt[ni] < EPSILONSTATE)//&&(!bdummy))// || (ni < numstates-1))) 
		       {
		          beta_nt[ni] = EPSILONSTATE;
		       }
		   }
		}
	
                ddenom = 0;
		
	        alpha_nt = alpha[nt];	     

		//gamma_nt - P(x=S| o_0,...,o_t)
	        //P(o_t+1,...,o_T|x_t=s,lambda) * P(o_0,...,o_t,xt=s|lambda)

		for (int ns = 0; ns < gamma_nt.length; ns++)
	        {
		   double dval = alpha_nt[ns]*beta_nt[ns];
	           ddenom += dval;
	           gamma_nt[ns] = dval;
		} 

		if (ddenom > 0)
		{
	           for (int ns = 0; ns < gamma_nt.length; ns++)
	           {
		       gamma_nt[ns] /= ddenom;
		   }
		}

                gammaObservedSum_combo_nt = gammaObservedSum[traindataObservedIndex[nt]];

                for (int ns = 0; ns < numstates; ns++)
		{
		   //first sum gamma over all common signatures
		   //updates probability of observing the signature when in the state
	           gammaObservedSum_combo_nt[ns] += gamma_nt[ns];
		}

		double dsum = 0;
		  
                //this compues the numerator portion
	        for (int ni = 0; ni < numstates; ni++)
	        {
		   double[] sumforsxi_ni = sumforsxi[ni]; //computing expected number of transition from state i
	           int[] transitionprobsindex_ni = transitionprobsindex[ni]; //indicies of non-zero transitions from state i
	           double[] transitionprobs_ni = transitionprobs[ni]; //probability of transitions from state i
	           int ntransitionprobsnum_ni = transitionprobsnum[ni]; //number of non-zero transitions from state i
	           double dalpha_nt_ni = alpha_nt[ni]; 
	           //sxi is P(q_t = S_i, q_(t+1) = S_j | O)

                   if (ntransitionprobsnum_ni < nsparsecutofflooser)
      	           {
		       //enough 0 transitionto use sparsity here
		       //looser cut off since the indirection is less of the total time
	              for (int nj = 0; nj < ntransitionprobsnum_ni; nj++)
		      {
		         int nmappedindex = transitionprobsindex_ni[nj];
		         //computes transition probability from state i to j
		         double dtempval = transitionprobs_ni[nmappedindex] *dalpha_nt_ni*tempproductbetaemiss[nmappedindex];
                         dsum += dtempval;
	      	         sumforsxi_ni[nmappedindex] = dtempval;
		      }
		   }
	           else
	           {
		      for (int nj = 0; nj < numstates; nj++)
		      {
		         //computes transition probability from state i to j
		         double dtempval = transitionprobs_ni[nj]*dalpha_nt_ni*tempproductbetaemiss[nj];
		         dsum += dtempval;
		         sumforsxi_ni[nj] = dtempval;
		      }
		   }
		}		 

		if (dsum > 0)
		{		   
		   //normalizing the numerator by the sum of the denominator and updating this iterations value for it
	           for (int ni = 0; ni < numstates; ni++)
	           {
		      int[] transitionprobsindex_ni = transitionprobsindex[ni];
	              double[] sumforsxi_ni = sumforsxi[ni];
	              double[] sxi_nseq_ni = sxi_nseq[ni];
	              int ntransitionprobsnum_ni = transitionprobsnum[ni];

                      if (ntransitionprobsnum_ni < nsparsecutoff)
	              {
		         //guessing sparse enough to avoid the indirections
		         for (int nj = 0; nj < ntransitionprobsnum_ni; nj++)
		         {
	      	            int nmappedindex = transitionprobsindex_ni[nj];
		            sxi_nseq_ni[nmappedindex] += sumforsxi_ni[nmappedindex]/dsum;
			 }
		      }
	              else
	              {
                         for (int nj = 0; nj < numstates; nj++)
	      	         {
			    sxi_nseq_ni[nj] += sumforsxi_ni[nj]/dsum;
		         }
		      }
		   }
	        }   	 
		beta_ntp1 = beta_nt;  //updating beta_ntp1 
	     }

	     double[] gammainitstore_nseq = gammainitstore[nseq];
	     for (int ns = 0; ns < numstates; ns++)
	     {
	        //storing the initial gamma from this iteration
	        gammainitstore_nseq[ns] = gamma_nt[ns];
	     }

	     //for (int nindex = 0; nindex < gammaObservedSum.length; nindex++)
	     for (int nindex = 0; nindex < nobserved; nindex++)
	     {
	        //going through all the gamma sufficient statistic

	        //only update for those combinations that were observed on this sequnce
	        //gets the observed combination and missing combination signatures
          	boolean[] traindataObservedValues_nindex = traindataObservedValues[nindex];
	        boolean[] traindataNotMissing_nindex = traindataNotMissing[nindex];

	           
	        double[] gammaObservedSum_nindex = gammaObservedSum[nindex];

	        for (int ns = 0; ns < numstates; ns++)
	        {
	           //going through each state
		   double[][] gammaksum_nseq_ns = gammaksum_nseq[ns];
	           double gammaObservedSum_nindex_ns = gammaObservedSum_nindex[ns];
	           for (int nmark = 0; nmark < numdatasets; nmark++)
	           {
		      //going through each mark
	      	      if (traindataNotMissing_nindex[nmark])
		      {
		         //only update non-missing
		         if (traindataObservedValues_nindex[nmark])
		         {
		            //updates the gamma sum for each mark when in state and observed 1
		            gammaksum_nseq_ns[nmark][1] += gammaObservedSum_nindex_ns;
		         }
		         else
		         {
		  	    //updates the gamma sum for each mark when in state and observed 0
		            gammaksum_nseq_ns[nmark][0] += gammaObservedSum_nindex_ns;
			 }
		      }
		   }
		}
	     }
	     

	  
	     //M step	      
	     if ((niteration >1) ||(nseq==(chromfiles.length-1)))
	     {
	        //executes the M-step after any pass through a sequence after one pass has been made through all sequences
	        double dsum = 0;

		//updating the inital probabilities
                for (int ni = 0; ni < numstates; ni++)
	        {
		   double dgammainitsum = 0;
                   for (int nitr = 0; nitr < chromfiles.length; nitr++)
	           {
	              dgammainitsum += gammainitstore[nitr][ni];
	           }

	           if ((bpseudo)&&(!bnopseudolast))
	           {
		      dgammainitsum++;
	           }

		   probinit[ni] = dgammainitsum;
                   dsum += dgammainitsum;
		}
 
	        for (int ni = 0; ni < numstates; ni++)
	        {
	           probinit[ni] /= dsum;		
	        }

		//this indicates if there is a change on the set of 0 probability transitions
	        boolean bchange = false;
	        for (int ni = 0; ni < transitionprobs.length; ni++)
                {
	           dsum = 0;
	           //computes the denominator for the transition probabilities

		   int[] transitionprobsindex_ni = transitionprobsindex[ni];
	           double[] transitionprobs_ni =  transitionprobs[ni];
		   int transitionprobsnum_ni = transitionprobsnum[ni];
	           for (int nj = 0; nj < transitionprobsnum_ni; nj++)
	           {
		      int ntransitionprobsindex_ni_nj = transitionprobsindex_ni[nj];
	              double dsxistoreitr = 0;
	              for (int nitr = 0; nitr < chromfiles.length; nitr++)
	              {
	                 dsxistoreitr += sxistore[nitr][ni][ntransitionprobsindex_ni_nj];
		      }

		      if (bpseudo)
		      {
		         //pseudo-count to avoid divide by zero
		         dsxistoreitr++;
		      }

		      transitionprobs_ni[ntransitionprobsindex_ni_nj] = dsxistoreitr;

	              dsum += dsxistoreitr;
		   }

	           for (int nj = 0; nj < transitionprobsnum_ni; nj++)
	           {
		      int ntransitionprobsindex_ni_nj = transitionprobsindex_ni[nj];
		      //computes the updated transition probabilities
		      transitionprobs_ni[ntransitionprobsindex_ni_nj] /= dsum;
	    
		      if ((transitionprobs_ni[ntransitionprobsindex_ni_nj] < dzerotransitioncutoff) && (ni != ntransitionprobsindex_ni_nj))
		      {
		         //if falls below threshold eliminate the transition probabilities
	       	         elim[ni][ntransitionprobsindex_ni_nj] = true;
		         bchange = true;
		         nelim++;
		         transitionprobs_ni[ntransitionprobsindex_ni_nj] = 0;
		      }
		   }
		}
	       
	        if (bchange)
	        {
	           //a transition was eliminated we need to update the probabilities
	           for (int ni = 0; ni < transitionprobs.length; ni++)
	           {
		      int nindex = 0;		        
                      ddenom = 0;
		      boolean[] elim_ni = elim[ni];
		      double[] transitionprobs_ni = transitionprobs[ni];
		      int[] transitionprobsindex_ni = transitionprobsindex[ni];
		      for (int nj = 0; nj < transitionprobs_ni.length; nj++)
	              {
		         if (!elim_ni[nj])
		         {
		            //we have not eliminated this transition
		            //stores its index in order and add sum to denominator
		            transitionprobsindex_ni[nindex] = nj;
		            ddenom += transitionprobs_ni[nj];
		            nindex++;
		         }
		      }

		      //renormalize the transition probabilities by the sum of the non-eliminated transitions
		      for (int nj = 0; nj < transitionprobs_ni.length; nj++)
	              {
		         transitionprobs_ni[nj] /= ddenom;
		      }
		      //update the number of valid transitions
		      transitionprobsnum[ni] = nindex; 
		   }

		   for (int ni = 0; ni < transitionprobs.length; ni++)
	           {
		      int nindex =0;
		      int[] transitionprobsindexCol_ni = transitionprobsindexCol[ni];
		      for (int nj = 0; nj < transitionprobs[ni].length; nj++)
	              {
		         if (!elim[nj][ni])
		         {
		            //copy into the column of i the index of all non-eliminated transitions of i
		            transitionprobsindexCol_ni[nindex] = nj;
		            nindex++;
		         }
		      }
		      //updates the number of non-zero transitions from column i
		      transitionprobsnumCol[ni] = nindex; 
		   }
		}
	    
	        //updating the emission parameters
	        for (int ns = 0; ns < numstates; ns++)
	        {
		   double[][] emissionprobs_ns = emissionprobs[ns];

	           for (int nmark = 0; nmark < emissionprobs_ns.length; nmark++)
	           {
		      double[] emissionprobs_ns_nmark = emissionprobs_ns[nmark];
		      //can't used a general gamma sum because of missing emission vals
		      double dgammadenom = 0;

	              //updates gamma sum
                      for (int nbucket = 0; nbucket < numbuckets; nbucket++)
		      {
		         emissionprobs_ns_nmark[nbucket] = 0;
	                
	                 for (int nitr = 0; nitr < chromfiles.length; nitr++)
	                 {
	                    emissionprobs_ns_nmark[nbucket] += gammaksumstore[nitr][ns][nmark][nbucket];
		         }

			 if ((bpseudo) && ((!bnopseudolast)||(nmark<emissionprobs_ns.length-1)))
			 {
			    emissionprobs_ns_nmark[nbucket]++;
			 }

		         dgammadenom += emissionprobs_ns_nmark[nbucket];
		      }


		      //added to avoid NA values
		      if (dgammadenom > 0)
		      {
		         for (int nbucket = 0; nbucket < numbuckets; nbucket++)
		         {
		            emissionprobs_ns_nmark[nbucket] /= dgammadenom;
		         }
		      }
		   }
		}
	     }

	     if (ChromHMM.BVERBOSE)
	     {
	        System.out.println("\t"+niteration+"\t"+dloglike);
	     }	     
	  }

	  double ddiff =(dloglike-dprevloglike);

	  //dconvergediff is only enforced if greater thanor equal to 0
          dprevloglike = dloglike;       

	  if (borderrows)
	  {
	     makeStateOrdering();
	  }

	  if (bordercols)
	  {
	     makeColOrdering();
	  }
	  //updates after each iteration the current status of the search
          printTransitionTable(niteration);
          printEmissionTable(niteration);
	  if (bprintimage)
	  {
             printEmissionImage(niteration);
             printTransitionImage(niteration);
	  }
	  printParametersToFile(niteration);

	  //we just completed a full iteration
          long ltimefinal =  System.currentTimeMillis();	  
	  double dtimechange = (ltimefinal-ltimeitr)/(double) 1000;
          bconverged = (((niteration >= nmaxiterations)||((ddiff< dconvergediff)&&(dconvergediff>=0)))||((dtimechange>nmaxseconds)&&(nmaxseconds>=0)));	  
          if (ChromHMM.BVERBOSE)
	  {
	     System.out.println(niteration+"\tTime Iteration\t"+dtimechange+"\t"+"\tElim\t"+nelim);
	     System.out.println("Full "+niteration+"\t"+dloglike+"\t"+dprevloglike+"\t"+ddiff);        
	  }


	  if (niteration == 1)
	  {
	      System.out.format("%10s %25s %10s %20s%n","Iteration","Estimated Log Likelihood", "Change","Total Time (secs)");
	      System.out.format("%10s %25s %10s %20s%n",""+niteration,""+nf3.format(dloglike),"-",""+nf1.format(dtimechange));
	  }
	  else
	  {
	      //System.out.format(niteration+"            "+nf3.format(dloglike)+"           "+nf3.format(ddiff)+"       "+nf1.format(dtimechange));
	      System.out.format("%10s %25s %10s %20s%n",""+niteration,""+nf3.format(dloglike),""+nf3.format(ddiff),""+nf1.format(dtimechange));
	  }
	  niteration++;
       }
       while (!bconverged);
    }

    //////////////////////////////////////////////////////////////////////////////////////////
    /**
     * This is the core procedure for learning the parameters of the model
     */
    public void trainParameters() throws IOException
    {
        NumberFormat nf3 = NumberFormat.getInstance();
        nf3.setMaximumFractionDigits(3);
	nf3.setGroupingUsed(false);
	nf3.setMinimumFractionDigits(3);

        NumberFormat nf1 = NumberFormat.getInstance();
	nf1.setMaximumFractionDigits(1);
	nf1.setMinimumFractionDigits(1);
	nf1.setGroupingUsed(false);

       int niteration = 1;

       boolean bconverged = false;

       double dzerotransitioncutoff = Math.pow(10,-nzerotransitionpower);

       //number of non-zero transition for the 
       int nsparsecutoff = (int) (numstates * ChromHMM.SPARSECUTOFFRATIO);

       //number of non-zero transition that need to be less than this at the looser cut-off
       int nsparsecutofflooser = (int) (numstates * ChromHMM.SPARSECUTOFFLOOSERRATIO);

       double dprevloglike;

       //stores the maximum number of locations in any sequence and in each sequence
       int[] numtime = new int[traindataObservedIndex.length];
       int nmaxtime = 0;
       for (int nseq = 0; nseq < traindataObservedIndex.length; nseq++)
       {
          numtime[nseq] = traindataObservedIndex[nseq].length;
          if (numtime[nseq] > nmaxtime)
	  {
      	     nmaxtime = numtime[nseq];
	  }
       }

       if (ChromHMM.BVERBOSE)
       {
          System.out.println("Maximum number of locations\t"+nmaxtime);
       }

       //stores the emission probability for the i^th combination of marks in the j^th state
       double[][] emissionproducts = new double[traindataObservedValues.length][numstates];

       //stores the scale for the emission probability for the i^th combination of marks
       double[] emissionproducts_scale = null;

       if (bscaleemissions)
       {
          emissionproducts_scale = new double[traindataObservedValues.length];
       }

       //stores temporary product terms
       double[] tempproductbetaemiss = new double[numstates];

       //This stores the alpha values at each time point and number of states
       double[][] alpha = new double[nmaxtime][numstates];

       //Temporary storage of the gamma's for each state
       double[] gamma_nt = new double[numstates];

       //Temporary storage of the beta values for each state
       double[] beta_nt = new double[numstates];

       //Temporary storage of the beta values for each state at the next time point
       double[] beta_ntp1 = new double[numstates];

       //stores the scaling value for each time point
       double[] scale = new double[nmaxtime];

       //stores the transition probabilities for each column
       double[][] coltransitionprobs = new double[numstates][numstates];

       //stores the sufficient statistic for the initital probability in each state for the last visit
       double[][] gammainitstore = new double[traindataObservedIndex.length][numstates];

       //stores the sufficient statistics for computing the transition probabilities cumulated for each iteration
       double[][][] sxistore = new double[traindataObservedIndex.length][numstates][numstates];

       //stores the sufficient statistic for computing the emission probabilities
       double[][][][] gammaksumstore = 
                 new double[traindataObservedIndex.length][numstates][numdatasets][numbuckets];

       //temporary storage in computation of sxi
       double[][] sumforsxi = new double[numstates][numstates];

       //stores the sum of the gamma values associated with each combination in each state
       double[][] gammaObservedSum = new double[traindataObservedValues.length][numstates];

       int nelim = 0;
       long ltimeitr= System.currentTimeMillis();
       dprevloglike = Double.NEGATIVE_INFINITY;

       do
       {
          dloglike = 0;	  

          for (int nseq = 0; nseq < traindataObservedIndex.length; nseq++)
          {
	      //going through each sequence

	     int[] traindataObservedIndex_nseq = traindataObservedIndex[nseq];
	     boolean[] traindataObservedSeqFlags_nseq = traindataObservedSeqFlags[nseq];

             double[][][] gammaksum_nseq = gammaksumstore[nseq];

	     for (int ns = 0; ns < gammaksum_nseq.length; ns++)
	     {
		 //resetting the gamma sufficient statistics in the current sequence
		double[][] gammaksum_nseq_ns = gammaksum_nseq[ns];
	        for (int nmark = 0; nmark < gammaksum_nseq_ns.length; nmark++)
	        {
		   for (int nbucket = 0; nbucket < numbuckets; nbucket++)
		   {
	              gammaksum_nseq_ns[nmark][nbucket] = 0;
	           }
	        }
	     }

	     double[][] sxi_nseq = sxistore[nseq];
	     for (int ni = 0; ni < sxi_nseq.length; ni++)
	     {
		 //reseeting the sxi sufficient statistics in the current sequence
		double[] sxi_nseq_ni = sxi_nseq[ni];
	        for (int nj = 0; nj < sxi_nseq_ni.length; nj++)
	        {
	           sxi_nseq_ni[nj] = 0;
	        }
	     }

	     //gammaObservedSum stores the weight for each combination of marks in each state
	     for (int ncombo = 0; ncombo < gammaObservedSum.length; ncombo++)
	     {
		 //resetting that to 0
		double[] gammaObservedSum_ncombo = gammaObservedSum[ncombo];
	        for (int ns = 0; ns < gammaObservedSum_ncombo.length; ns++)
	        {
	           gammaObservedSum_ncombo[ns] = 0;
	        }
	     }

             if (bscaleemissions)
	     {
	        for (int ni = 0; ni < emissionproducts.length; ni++)
		{
		   //going through each combination of marks
		   if (traindataObservedSeqFlags_nseq[ni])
		   {
		      //this signature of marks is observed on the current chromosome so
		      //updating its emission probabilities
		       double[] emissionproducts_ni = emissionproducts[ni];
		       boolean[] traindataObservedValues_ni = traindataObservedValues[ni];
		       boolean[] traindataNotMissing_ni = traindataNotMissing[ni];

		       for (int ns = 0; ns < numstates; ns++)
		       {
	       	          emissionproducts_ni[ns] = 1;
		       }
		       emissionproducts_scale[ni] = 0;

		       for (int nmod = 0; nmod < numdatasets; nmod++)
		       {
		          for (int ns = 0; ns < numstates; ns++)
			  {
			     if (traindataNotMissing_ni[nmod])
			     {
			        //we are include this marks emission probability
				if (traindataObservedValues_ni[nmod])
			        {
				   emissionproducts_ni[ns] *= emissionprobs[ns][nmod][1];
				}
				else
				{
				    emissionproducts_ni[ns] *= emissionprobs[ns][nmod][0];
				}
			     }
			     // otherwise treated as missing omitting from product
			  }

			  double dmaxval = 0;
			  for (int ns = 0; ns < numstates; ns++)
			  {
			     if (emissionproducts_ni[ns] > dmaxval)
			     {
			        dmaxval = emissionproducts_ni[ns];
			     }
			  }

			  if (dmaxval <= 0)//EPSILONEMISSIONS)
			  {
			     for (int ns = 0; ns < numstates; ns++)
			     {
			        emissionproducts_ni[ns] = 1;
			     }
			     //approximating log-likelihood
			     emissionproducts_scale[ni] += Math.log(EPSILONEMISSIONS);
			  }
			  else
			  {
			     for (int ns = 0; ns < numstates; ns++)
			     {
			        emissionproducts_ni[ns]/= dmaxval;
			     }
			     emissionproducts_scale[ni] += Math.log(dmaxval);
			  }
		       }
		   }
		}
	     }
	     else
	     {
	        for (int ni = 0; ni < emissionproducts.length; ni++)
	        {
	           //going through each combination of marks
	           if (traindataObservedSeqFlags_nseq[ni])
		   {
		      //this signature of marks is observed on the current chromosome so
		      //updating its emission probabilities
		      double[] emissionproducts_ni = emissionproducts[ni];
		      boolean[] traindataObservedValues_ni = traindataObservedValues[ni];
		      boolean[] traindataNotMissing_ni = traindataNotMissing[ni];		  

		      boolean ballzero = true;

	              for (int ns = 0; ns < numstates; ns++)
	              {
		         double dproduct = 1;
		         double[][] emissionprobs_ns = emissionprobs[ns];

		         for (int nmod = 0; nmod < numdatasets; nmod++)
		         {
			    if (traindataNotMissing_ni[nmod])
			    {
			       //we are include this marks emission probability
		               if (traindataObservedValues_ni[nmod])
		               {
				  //System.out.println("positive\t"+ns+"\t"+nmod+"\t1\t"+emissionprobs_ns[nmod][1]);
		                  dproduct *= emissionprobs_ns[nmod][1];
			       }
		               else 
		               {
				  ///System.out.println("negative\t"+ns+"\t"+nmod+"\t0\t"+emissionprobs_ns[nmod][0]);
		                  dproduct *= emissionprobs_ns[nmod][0];
			       }
			    }
			    // otherwise treated as missing omitting from product
			 }
		         //System.out.println(ns+"\t"+dproduct);
		         emissionproducts_ni[ns] = dproduct;

		         if (dproduct >= EPSILONEMISSIONS)
		         {
		            ballzero = false;
			 }
		      }

		      if (ballzero)
		      {
		         for (int ns = 0; ns < numstates; ns++)
		         {
			    emissionproducts_ni[ns] = EPSILONEMISSIONS;
			 }
		      }
		   }
		}
	     }

	     //initial probability in state s is initial probability times emission probability at first position
	     double[] alpha_nt = alpha[0];
	     double dscale = 0;
	     double[] emissionproducts_nobserveindex =emissionproducts[traindataObservedIndex_nseq[0]];

 	     for (int ns = 0; ns < numstates; ns++)
             {
	        alpha_nt[ns] = probinit[ns] * emissionproducts_nobserveindex[ns];
	        dscale += alpha_nt[ns];
	     }
	     
	     scale[0] = dscale;

	     //alpha_t(s)=P(o_0,...,o_t,x_t=s|lambda)
	     //converts the alpha terms to probabilities
	     if (bscalebeta)
	     {
	        if (dscale == 0)
	        {
		   for (int ns = 0; ns < numstates; ns++)
		   {
		      //added
		      if ((alpha_nt[ns] < EPSILONSTATE)&&(emissionproducts_nobserveindex[ns]>0))//(ns < numstates-1))) 
		      {
		         alpha_nt[ns] = EPSILONSTATE;
	       	         dscale += EPSILONSTATE;
		      }
		   }
	           scale[0] = dscale;
		}

	        for (int ns = 0; ns < numstates; ns++)
	        {
                   alpha_nt[ns] /= dscale;

		   if ((alpha_nt[ns] < EPSILONSTATE)&&(emissionproducts_nobserveindex[ns]>0))// (ns < numstates-1))) 
		   {
		      alpha_nt[ns] = EPSILONSTATE;
		   }
		}
	     }
      	     else
	     {
	        for (int ns = 0; ns < numstates; ns++)
	        {
                   alpha_nt[ns] /= dscale;
		}
	     }


	     dloglike += Math.log(dscale); 

             if (bscaleemissions)
	     {
	        dloglike += emissionproducts_scale[traindataObservedIndex_nseq[0]];
       	     }

	     //stores in coltransitionprobs the transpose of transitionprobs
	     for (int ni = 0; ni < numstates; ni++)
	     {
	        double[] coltransitionprobs_ni = coltransitionprobs[ni];
	        for (int nj = 0; nj < numstates; nj++)
	        {
	           coltransitionprobs_ni[nj] = transitionprobs[nj][ni];
	        }
	     }

	     //forward step
	     int numtime_nseq = numtime[nseq];
	     for (int nt = 1; nt < numtime_nseq; nt++)
	     {
	        //the actual observed combination at position t	        
	        double[] alpha_ntm1 = alpha[nt-1];
	        alpha_nt = alpha[nt];
	      
	        dscale = 0;
		emissionproducts_nobserveindex = emissionproducts[traindataObservedIndex_nseq[nt]];
	        for (int ns = 0; ns < numstates; ns++)
	        {
		   //stores the emission product for each location on the chromosome		   

		   int transitionprobsnumCol_ns = transitionprobsnumCol[ns];
		   int[] transitionprobsindexCol_ns = transitionprobsindexCol[ns];
	           double[] coltransitionprobs_ns = coltransitionprobs[ns];

	           double dtempsum = 0;
                   if (transitionprobsnumCol_ns < nsparsecutoff)
		   {
		       //number of transitions is sparse enough worth going through the extra redirection
	              for (int nj = 0; nj < transitionprobsnumCol_ns; nj++)
		      {
		         //for each next state computing inner sum of all previous alpha and the transition probability
		         //for all non-zero transitions into the state
		         int nmappedindex = transitionprobsindexCol_ns[nj];
		         dtempsum += coltransitionprobs_ns[nmappedindex]*alpha_ntm1[nmappedindex];
		      }
		   }
		   else
	           {
		       //avoid the redirect and multiply by 0
                      for (int nj = 0; nj < numstates; nj++)
		      {
		         dtempsum += coltransitionprobs_ns[nj]*alpha_ntm1[nj];
		      }
		   }
	           //multiply the transition sum by the emission probability
		   double dalphaval = dtempsum*emissionproducts_nobserveindex[ns];
	           alpha_nt[ns] = dalphaval;
		   //System.out.println(ns+"\t"+alpha_nt[ns]+"\t"+dtempsum+"\t"+emissionproducts_nobserveindex[ns]);
	           dscale += dalphaval;
	        }


		//rescaling alpha
	        scale[nt] = dscale;
                //scale_t(s)=P(o_0,...,o_t|lambda) summed over all states

		if (bscalebeta)
		{
	           if (dscale == 0)
	           {
		      for (int ns = 0; ns < numstates; ns++)
		      {
		         //added
		         if ((alpha_nt[ns] < EPSILONSTATE)&&(emissionproducts_nobserveindex[ns]>0))//(ns < numstates-1))) 
		         {
		            alpha_nt[ns] = EPSILONSTATE;
			    dscale += EPSILONSTATE;
		         }
		      }
		      scale[nt] = dscale;
		   }

	           for (int ns = 0; ns < numstates; ns++)
	           {
                      alpha_nt[ns] /= dscale;

		      //if ((alpha_nt[ns] < EPSILONSTATE)&&(!bdummy || (ns < numstates-1))) 
		      if ((alpha_nt[ns] < EPSILONSTATE)&&(emissionproducts_nobserveindex[ns]>0))
		      {
		         alpha_nt[ns] = EPSILONSTATE;
		      }
		   }
		}
		else
		{
	           for (int ns = 0; ns < numstates; ns++)
	           {
                      alpha_nt[ns] /= dscale;
		   }
		}

      	        dloglike += Math.log(dscale);

                if (bscaleemissions)
		{
		   dloglike += emissionproducts_scale[traindataObservedIndex_nseq[nt]];
	        }

	     }
	   

	     //backward step
	     //beta_t(s)=P(o_t+1,...,o_T|x_t=s,lambda)
	     int nlastindex = numtime_nseq-1;
	     double dinitval;
	     if (bscalebeta)
	     {
		dinitval = 1.0/numstates;
	     }
	     else
	     {
                dinitval = 1.0/scale[nlastindex];
	     }

             for (int ns = 0; ns < numstates; ns++)
	     {
	        beta_ntp1[ns] = dinitval;
	     }

	     double ddenom = 0;	       

	     //gamma_nt - P(x=S| o_0,...,o_t)
	     //P(o_t+1,...,o_T|x_t=s,lambda) * P(o_0,...,o_t,xt=s|lambda)
	     alpha_nt = alpha[nlastindex];
	     for (int ns = 0; ns < gamma_nt.length; ns++)
	     {
	        double dval = alpha_nt[ns]*beta_ntp1[ns];
	        ddenom += dval;
	        gamma_nt[ns] = dval;
	     }

	     if (ddenom > 0)
	     {
                for (int ns = 0; ns < gamma_nt.length; ns++)
	        {
		   gamma_nt[ns] /=ddenom;
		}
	     }

	     double[] gammaObservedSum_combo_nt = gammaObservedSum[traindataObservedIndex_nseq[nlastindex]];
		
	     for (int ns = 0; ns < numstates; ns++)
	     { 
                //first sum gamma over all common signatures		     
	        //updates probability of observing the signature when in the state
	        gammaObservedSum_combo_nt[ns] += gamma_nt[ns];
	     }

	     for (int nt = nlastindex - 1; nt >= 0; nt--)
	     {
	        int ntp1 = (nt+1);
		   
	        double[] emissionproducts_combo_ntp1 = emissionproducts[traindataObservedIndex_nseq[ntp1]];

		for (int ns = 0; ns < numstates; ns++)
	        {
	           tempproductbetaemiss[ns] = beta_ntp1[ns]*emissionproducts_combo_ntp1[ns];		   
		}

	        if (bscaleemissions)
	        {
		   //adding here to help numerical stability
		   double dmaxval = 0;
		   for (int ns = 0; ns < numstates; ns++)
		   {
		      if (tempproductbetaemiss[ns] > dmaxval)
		      {
			  dmaxval = tempproductbetaemiss[ns];
		      }
		   }
                   for (int ns = 0; ns < numstates; ns++)
		   {
		      tempproductbetaemiss[ns] = tempproductbetaemiss[ns]/dmaxval;
		   }
		}

		//double dscaleinv = 1.0/scale[nt];
		double dsumbeta = 0;
		double dscale_nt = scale[nt];
                //scale_t(s)=P(o_0,...,o_t|lambda) summed over all states

	        for (int ni = 0; ni < numstates; ni++)
		{
	           double dtempsum = 0;
	           int[] transitionprobsindex_ni =  transitionprobsindex[ni];
	           double[] transitionprobs_ni = transitionprobs[ni];
		   int transitionprobsnum_ni = transitionprobsnum[ni];

                   if (transitionprobsnum_ni < nsparsecutoff)
		   {
		       //sparse enought to pay the indirection penalty
		      for (int nj = 0; nj < transitionprobsnum_ni; nj++)
		      {
		         //for each state summing over transition probability to state j, emission probablity in j at next step
		         //and probability of observing the remaining sequence
			 int nmappedindexouter = transitionprobsindex_ni[nj];
		         dtempsum += transitionprobs_ni[nmappedindexouter]*tempproductbetaemiss[nmappedindexouter];
		      }
		   }
		   else
		   {
		       //not trying to exploit sparsity here
                      for (int nj = 0; nj < numstates; nj++)
		      {
		         //for each state summing over transition probability to state j, emission probablity in j at next step
		         //and probability of observing the remaining sequence
			 dtempsum += transitionprobs_ni[nj]*tempproductbetaemiss[nj];
		      }
		   }



		   if (bscalebeta)
		   {
		       beta_nt[ni] = dtempsum;
		       dsumbeta += dtempsum;
		   }
		   else
		   {
		      double dratio = dtempsum/dscale_nt;
		      if (dratio > Double.MAX_VALUE)
		      {
		         beta_nt[ni] = Double.MAX_VALUE;//dtempsum/dscale_nt;
		      }
		      else
		      {
		         beta_nt[ni] = dratio;
		      }
		   }
		}

		
		if (bscalebeta)
		{
                   for (int ni = 0; ni < numstates; ni++)
		   {
		      beta_nt[ni]/= dsumbeta;

		      if (beta_nt[ni] < EPSILONSTATE)//&& (!bdummy))// || !ballemisszero))// (!bdummy || (ni < numstates-1)))
		      { 
		         beta_nt[ni] = EPSILONSTATE;
		      }
		   }
		}
				
	     
                ddenom = 0;
		
	        alpha_nt = alpha[nt];	     

		//gamma_nt - P(x=S| o_0,...,o_t)
	        //P(o_t+1,...,o_T|x_t=s,lambda) * P(o_0,...,o_t,xt=s|lambda)

		for (int ns = 0; ns < gamma_nt.length; ns++)
	        {
		   double dval = alpha_nt[ns]*beta_nt[ns];
	           ddenom += dval;
	           gamma_nt[ns] = dval;
		} 

		if (ddenom > 0)
		{
	           for (int ns = 0; ns < gamma_nt.length; ns++)
	           {
		      gamma_nt[ns] /= ddenom;
		   }
		}

                gammaObservedSum_combo_nt = gammaObservedSum[traindataObservedIndex_nseq[nt]];

                for (int ns = 0; ns < numstates; ns++)
		{
		   //first sum gamma over all common signatures
		   //updates probability of observing the signature when in the state
	           gammaObservedSum_combo_nt[ns] += gamma_nt[ns];
		}

		double dsum = 0;
		  
                //this compues the numerator portion
	        for (int ni = 0; ni < numstates; ni++)
	        {
		   double[] sumforsxi_ni = sumforsxi[ni]; //computing expected number of transition from state i
	           int[] transitionprobsindex_ni = transitionprobsindex[ni]; //indicies of non-zero transitions from state i
	           double[] transitionprobs_ni = transitionprobs[ni]; //probability of transitions from state i
	           int ntransitionprobsnum_ni = transitionprobsnum[ni]; //number of non-zero transitions from state i
	           double dalpha_nt_ni = alpha_nt[ni]; 
	           //sxi is P(q_t = S_i, q_(t+1) = S_j | O)

                   if (ntransitionprobsnum_ni < nsparsecutofflooser)
      	           {
		       //enough 0 transitionto use sparsity here
		       //looser cut off since the indirection is less of the total time
	              for (int nj = 0; nj < ntransitionprobsnum_ni; nj++)
		      {
		         int nmappedindex = transitionprobsindex_ni[nj];
		         //computes transition probability from state i to j
		         double dtempval = transitionprobs_ni[nmappedindex] *dalpha_nt_ni*tempproductbetaemiss[nmappedindex];
                         dsum += dtempval;
	      	         sumforsxi_ni[nmappedindex] = dtempval;
		      }
		   }
	           else
	           {
		      for (int nj = 0; nj < numstates; nj++)
		      {
		         //computes transition probability from state i to j
		         double dtempval = transitionprobs_ni[nj]*dalpha_nt_ni*tempproductbetaemiss[nj];
		         dsum += dtempval;
		         sumforsxi_ni[nj] = dtempval;
		       }
		   }
		}		 
		   
		if (dsum > 0)
		{
		   //normalizing the numerator by the sum of the denominator and updating this iterations value for it
	           for (int ni = 0; ni < numstates; ni++)
	           {
		      int[] transitionprobsindex_ni = transitionprobsindex[ni];
	              double[] sumforsxi_ni = sumforsxi[ni];
	              double[] sxi_nseq_ni = sxi_nseq[ni];
	              int ntransitionprobsnum_ni = transitionprobsnum[ni];

                      if (ntransitionprobsnum_ni < nsparsecutoff)
	              {
		         //guessing sparse enough to avoid the indirections
		         for (int nj = 0; nj < ntransitionprobsnum_ni; nj++)
		         {
	      	            int nmappedindex = transitionprobsindex_ni[nj];
		            sxi_nseq_ni[nmappedindex] += sumforsxi_ni[nmappedindex]/dsum;
			 }
		      }
	              else
	              {
                         for (int nj = 0; nj < numstates; nj++)
	      	         {
			    sxi_nseq_ni[nj] += sumforsxi_ni[nj]/dsum;
			 }
		      }
		   }
	        }   	 
		beta_ntp1 = beta_nt;  //updating beta_ntp1 
	     }

	     double[] gammainitstore_nseq = gammainitstore[nseq];
	     for (int ns = 0; ns < numstates; ns++)
	     {
	        //storing the initial gamma from this iteration
	        gammainitstore_nseq[ns] = gamma_nt[ns];
	     }

	     for (int nindex = 0; nindex < gammaObservedSum.length; nindex++)
	     {
		 //going through all the gamma sufficient statistic
	        if (traindataObservedSeqFlags_nseq[nindex])
	        {
	           //only update for those combinations that were observed on this sequnce
		   //gets the observed combination and missing combination signatures
		   boolean[] traindataObservedValues_nindex = traindataObservedValues[nindex];
		   boolean[] traindataNotMissing_nindex = traindataNotMissing[nindex];

	           
	           double[] gammaObservedSum_nindex = gammaObservedSum[nindex];

		   for (int ns = 0; ns < numstates; ns++)
	           {
		      //going through each state
		      double[][] gammaksum_nseq_ns = gammaksum_nseq[ns];
	              double gammaObservedSum_nindex_ns = gammaObservedSum_nindex[ns];
		      for (int nmark = 0; nmark < numdatasets; nmark++)
	              {
			  //going through each mark
			  if (traindataNotMissing_nindex[nmark])
		         {
			     //only update non-missing
		            if (traindataObservedValues_nindex[nmark])
		            {
		      	       //updates the gamma sum for each mark when in state and observed 1
		               gammaksum_nseq_ns[nmark][1] += gammaObservedSum_nindex_ns;
			    }
		            else
		            {
		       	       //updates the gamma sum for each mark when in state and observed 0
		               gammaksum_nseq_ns[nmark][0] += gammaObservedSum_nindex_ns;
			    }
			 }
		      }
		   }
		}
	     }

	  
	     //M step	      
	     if ((niteration >1) ||(nseq==(traindataObservedIndex.length-1)))
	     {
	        //executes the M-step after any pass through a sequence after one pass has been made through all sequences
	        double dsum = 0;

		//updating the inital probabilities
                for (int ni = 0; ni < numstates; ni++)
	        {
		   double dgammainitsum = 0;
                   for (int nitr = 0; nitr < traindataObservedIndex.length; nitr++)
	           {
	              dgammainitsum += gammainitstore[nitr][ni];
	           }

	           if ((bpseudo)&&(!bnopseudolast))
	           {
		      dgammainitsum++;
		   }

		   probinit[ni] = dgammainitsum;
                   dsum += dgammainitsum;
		}
 
	        for (int ni = 0; ni < numstates; ni++)
	        {
	           probinit[ni] /= dsum;		
	        }

		//this indicates if there is a change on the set of 0 probability transitions
	        boolean bchange = false;
	        for (int ni = 0; ni < transitionprobs.length; ni++)
                {
	           dsum = 0;
	           //computes the denominator for the transition probabilities

		   int[] transitionprobsindex_ni = transitionprobsindex[ni];
	           double[] transitionprobs_ni =  transitionprobs[ni];
		   int transitionprobsnum_ni = transitionprobsnum[ni];
	           for (int nj = 0; nj < transitionprobsnum_ni; nj++)
	           {
		      int ntransitionprobsindex_ni_nj = transitionprobsindex_ni[nj];
	              double dsxistoreitr = 0;
	              for (int nitr = 0; nitr < traindataObservedIndex.length; nitr++)
	              {
	                 dsxistoreitr += sxistore[nitr][ni][ntransitionprobsindex_ni_nj];
		      }

		      if (bpseudo)
		      {
		         //pseudo-count to avoid divide by zero
		         dsxistoreitr++;
		      }

		      transitionprobs_ni[ntransitionprobsindex_ni_nj] = dsxistoreitr;

	              dsum += dsxistoreitr;
		   }

	           for (int nj = 0; nj < transitionprobsnum_ni; nj++)
	           {
		      int ntransitionprobsindex_ni_nj = transitionprobsindex_ni[nj];
		      //computes the updated transition probabilities
		      transitionprobs_ni[ntransitionprobsindex_ni_nj] /= dsum;
	    
		      if ((transitionprobs_ni[ntransitionprobsindex_ni_nj] < dzerotransitioncutoff) && (ni != ntransitionprobsindex_ni_nj))
		      {
		         //if falls below threshold eliminate the transition probabilities
	       	         elim[ni][ntransitionprobsindex_ni_nj] = true;
		         bchange = true;
		         nelim++;
		         transitionprobs_ni[ntransitionprobsindex_ni_nj] = 0;
		      }
		   }
		}
	       
	        if (bchange)
	        {
	           //a transition was eliminated we need to update the probabilities
	           for (int ni = 0; ni < transitionprobs.length; ni++)
	           {
		      int nindex = 0;		        
                      ddenom = 0;
		      boolean[] elim_ni = elim[ni];
		      double[] transitionprobs_ni = transitionprobs[ni];
		      int[] transitionprobsindex_ni = transitionprobsindex[ni];
		      for (int nj = 0; nj < transitionprobs_ni.length; nj++)
	              {
		         if (!elim_ni[nj])
		         {
		            //we have not eliminated this transition
		            //stores its index in order and add sum to denominator
		            transitionprobsindex_ni[nindex] = nj;
		            ddenom += transitionprobs_ni[nj];
		            nindex++;
		         }
		      }

		      //renormalize the transition probabilities by the sum of the non-eliminated transitions
		      for (int nj = 0; nj < transitionprobs_ni.length; nj++)
	              {
		         transitionprobs_ni[nj] /= ddenom;
		      }
		      //update the number of valid transitions
		      transitionprobsnum[ni] = nindex; 
		   }

		   for (int ni = 0; ni < transitionprobs.length; ni++)
	           {
		      int nindex =0;
		      int[] transitionprobsindexCol_ni = transitionprobsindexCol[ni];
		      for (int nj = 0; nj < transitionprobs[ni].length; nj++)
	              {
		         if (!elim[nj][ni])
		         {
		            //copy into the column of i the index of all non-eliminated transitions of i
		            transitionprobsindexCol_ni[nindex] = nj;
		            nindex++;
		         }
		      }
		      //updates the number of non-zero transitions from column i
		      transitionprobsnumCol[ni] = nindex; 
		   }
		}
	    
	        //updating the emission parameters
	        for (int ns = 0; ns < numstates; ns++)
	        {
		   double[][] emissionprobs_ns = emissionprobs[ns];

	           for (int nmark = 0; nmark < emissionprobs_ns.length; nmark++) 
	           {
		      double[] emissionprobs_ns_nmark = emissionprobs_ns[nmark];
		      //can't used a general gamma sum because of missing emission vals
		      double dgammadenom = 0;

	              //updates gamma sum
                      for (int nbucket = 0; nbucket < numbuckets; nbucket++)
		      {
		         emissionprobs_ns_nmark[nbucket] = 0;
	                
	                 for (int nitr = 0; nitr < traindataObservedIndex.length; nitr++)
	                 {
	                    emissionprobs_ns_nmark[nbucket] += gammaksumstore[nitr][ns][nmark][nbucket];
		         }


			 if ((bpseudo) && ((!bnopseudolast)||(nmark<emissionprobs_ns.length-1)))
		         {
			    emissionprobs_ns_nmark[nbucket]++;
		         }
		         dgammadenom += emissionprobs_ns_nmark[nbucket];
		      }

		      //added to avoid NA values
		      if (dgammadenom > 0)
		      {
		         for (int nbucket = 0; nbucket < numbuckets; nbucket++)
		         {
		            emissionprobs_ns_nmark[nbucket] /= dgammadenom;
			 }
		      }
		   }
		}
	     }

	     if (ChromHMM.BVERBOSE)
	     {
	        System.out.println("\t"+niteration+"\t"+dloglike);
	     }	     
	  }

	  double ddiff =(dloglike-dprevloglike);

	  //dconvergediff is only enforced if greater thanor equal to 0
          dprevloglike = dloglike;       

	  if (borderrows)
	  {
	     makeStateOrdering();
	  }

	  if (bordercols)
	  {
	     makeColOrdering();
	  }
	  //updates after each iteration the current status of the search
          printTransitionTable(niteration);
          printEmissionTable(niteration);
	  if (bprintimage)
	  {
             printEmissionImage(niteration);
             printTransitionImage(niteration);
	  }
	  printParametersToFile(niteration);

	  //we just completed a full iteration
          long ltimefinal =  System.currentTimeMillis();	  
	  double dtimechange = (ltimefinal-ltimeitr)/(double) 1000;
          bconverged = (((niteration >= nmaxiterations)||((ddiff< dconvergediff)&&(dconvergediff>=0)))||((dtimechange>nmaxseconds)&&(nmaxseconds>=0)));	  
          if (ChromHMM.BVERBOSE)
	  {
	     System.out.println(niteration+"\tTime Iteration\t"+dtimechange+"\t"+"\tElim\t"+nelim);
	     System.out.println("Full "+niteration+"\t"+dloglike+"\t"+dprevloglike+"\t"+ddiff);        
	  }


	  if (niteration == 1)
	  {
	      System.out.format("%10s %25s %10s %20s%n","Iteration","Estimated Log Likelihood", "Change","Total Time (secs)");
	      System.out.format("%10s %25s %10s %20s%n",""+niteration,""+nf3.format(dloglike),"-",""+nf1.format(dtimechange));
	  }
	  else
	  {
	      //System.out.format(niteration+"            "+nf3.format(dloglike)+"           "+nf3.format(ddiff)+"       "+nf1.format(dtimechange));
	      System.out.format("%10s %25s %10s %20s%n",""+niteration,""+nf3.format(dloglike),""+nf3.format(ddiff),""+nf1.format(dtimechange));
	  }
	  niteration++;
       }
       while (!bconverged);
    }


    /////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create a new thread.
    class NewThread implements Runnable 
    {
       
	int[] traindataObservedIndex_nseq;
	boolean[] traindataObservedSeqFlags_nseq;
	double[][][] gammaksum_nseq;
	double[][] sxi_nseq;
	int numtime_nseq;
	double[] gammainitstore_nseq;
	double[][][] gammaObservedSum_Pool;
	double[][][] alpha_Pool;
	double[][] emissionproducts;
	double[][] gamma_nt_Pool;
	double[][] coltransitionprobs;
	double[][] scale_Pool;
	double[][] beta_nt_Pool;
	double[][] beta_ntp1_Pool;
	double[][] tempproductbetaemiss_Pool;
	int nsparsecutoff;
	int nsparsecutofflooser;
	double[][][] sumforsxi_Pool;
	double[] dloglikeA;
	int nseq;
        double[] emissionproducts_scale;

	NewThread(int[] traindataObservedIndex_nseq,
		  boolean[] traindataObservedSeqFlags_nseq,
		  double[][][] gammaksum_nseq,
		  double[][] sxi_nseq,
		  int numtime_nseq,
		  double[] gammainitstore_nseq,
		  double[][][] gammaObservedSum_Pool,
		  double[][][] alpha_Pool,
		  double[][] emissionproducts,
		  double[][] gamma_nt_Pool,
		  double[][] coltransitionprobs,
		  double[][] scale_Pool,
		  double[][] beta_nt_Pool,
		  double[][] beta_ntp1_Pool,
		  double[][] tempproductbetaemiss_Pool,
		  double[][][] sumforsxi_Pool,
		  int nsparsecutoff,
		  int nsparsecutofflooser,
		  double[] dloglikeA,
		  int nseq,
		  double[] emissionproducts_scale)
        {
	    this.traindataObservedIndex_nseq = traindataObservedIndex_nseq;
	    this.traindataObservedSeqFlags_nseq = traindataObservedSeqFlags_nseq;
	    this.gammaksum_nseq = gammaksum_nseq;
	    this.sxi_nseq = sxi_nseq;
	    this.numtime_nseq = numtime_nseq;
	    this.gammainitstore_nseq = gammainitstore_nseq;
	    this.gammaObservedSum_Pool = gammaObservedSum_Pool;
	    this.alpha_Pool = alpha_Pool;
	    this.emissionproducts = emissionproducts;
	    this.gamma_nt_Pool = gamma_nt_Pool;
	    this.coltransitionprobs = coltransitionprobs;
	    this.scale_Pool = scale_Pool;
	    this.beta_nt_Pool = beta_nt_Pool;
	    this.beta_ntp1_Pool = beta_ntp1_Pool;
	    this.tempproductbetaemiss_Pool = tempproductbetaemiss_Pool;
	    this.nsparsecutoff = nsparsecutoff;
	    this.nsparsecutofflooser = nsparsecutofflooser;
	    this.sumforsxi_Pool = sumforsxi_Pool;
	    this.dloglikeA = dloglikeA;
	    this.nseq = nseq;
	    this.emissionproducts_scale = emissionproducts_scale;
	}

	private int slotavailable()
	{
	   for (int ni = 0; ni < threadslots.length; ni++)
	   {
	      if (!threadslots[ni])
	      {
	         threadslots[ni] = true;
	         return ni;
	      }
	   }
	   return -1;	    
	}

	// This is the entry point for the second thread.
	public void run() 
        {
	    int nprocess;

	    synchronized (objlock)
	    {
	       while ((nprocess=slotavailable())==-1)
	       {
                  try
                  {
	             objlock.wait();
	          }
	          catch (InterruptedException ex) {}
	       }
	    }

	    estep(gammaObservedSum_Pool[nprocess],
		  alpha_Pool[nprocess],
		  gamma_nt_Pool[nprocess],
		  scale_Pool[nprocess],
		  beta_nt_Pool[nprocess],
		  beta_ntp1_Pool[nprocess],
		  tempproductbetaemiss_Pool[nprocess],
		  sumforsxi_Pool[nprocess]);
	    
	    synchronized(objlock)
	    {
		nlaunched--;
		threadslots[nprocess] = false;
		objlock.notifyAll();
	    }

	}

        public void estep(
		  double[][] gammaObservedSum,
		  double[][] alpha,
		  double[] gamma_nt,
		  double[] scale,
		  double[] beta_nt,
		  double[] beta_ntp1,
		  double[] tempproductbetaemiss,
		  double[][] sumforsxi)
	{ 
	      //going through each sequence

	   double dloglikeseq = 0;

	   for (int ns = 0; ns < gammaksum_nseq.length; ns++)
	   {
	      //resetting the gamma sufficient statistics in the current sequence
	      double[][] gammaksum_nseq_ns = gammaksum_nseq[ns];
	      for (int nmark = 0; nmark < gammaksum_nseq_ns.length; nmark++)
	      {
	         for (int nbucket = 0; nbucket < numbuckets; nbucket++)
	         {
	            gammaksum_nseq_ns[nmark][nbucket] = 0;
	         }
	      }
           }

	   for (int ni = 0; ni < sxi_nseq.length; ni++)
	   {
	      //reseeting the sxi sufficient statistics in the current sequence
	      double[] sxi_nseq_ni = sxi_nseq[ni];
	      for (int nj = 0; nj < sxi_nseq_ni.length; nj++)
	      {
                 sxi_nseq_ni[nj] = 0;
	      }
	   }

	   //gammaObservedSum stores the weight for each combination of marks in each state
	   for (int ncombo = 0; ncombo < gammaObservedSum.length; ncombo++)
	   {
	      //resetting that to 0
	      double[] gammaObservedSum_ncombo = gammaObservedSum[ncombo];
	      for (int ns = 0; ns < gammaObservedSum_ncombo.length; ns++)
	      {
	         gammaObservedSum_ncombo[ns] = 0;
	      }
	   }

	   //initial probability in state s is initial probability times emission probability at first position
	   double[] alpha_nt = alpha[0];
           double dscale = 0;
	   double[] emissionproducts_nobserveindex =emissionproducts[traindataObservedIndex_nseq[0]];
 	   for (int ns = 0; ns < numstates; ns++)
           {
	      alpha_nt[ns] = probinit[ns] * emissionproducts_nobserveindex[ns];
	      dscale += alpha_nt[ns];
	      
	   }
	   scale[0] = dscale;

	   //alpha_t(s)=P(o_0,...,o_t,x_t=s|lambda)
	   //converts the alpha terms to probabilities

	   if (bscalebeta)
	   {
	      if (dscale == 0)
	      {
		 for (int ns = 0; ns < numstates; ns++)
		 {
		    //added
		    if ((alpha_nt[ns] < EPSILONSTATE)&&(emissionproducts_nobserveindex[ns]>0))//(ns < numstates-1))) 
		    {
		       alpha_nt[ns] = EPSILONSTATE;
		       dscale += EPSILONSTATE;
	            }
       	         }
	         scale[0] = dscale;
	      }

	      for (int ns = 0; ns < numstates; ns++)
	      {
                 alpha_nt[ns] /= dscale;

	         if ((alpha_nt[ns] < EPSILONSTATE)&&(emissionproducts_nobserveindex[ns]>0))//(ns < numstates-1))) 
       	         {
	            alpha_nt[ns] = EPSILONSTATE;
	         }
	      }
	   }
      	   else
	   {
	      for (int ns = 0; ns < numstates; ns++)
	      {
                 alpha_nt[ns] /= dscale;
	      }
	   }
	   //for (int ni = 0; ni < numstates; ni++)
	   //{
	   //   alpha_nt[ni] /= dscale;
	   //}

	   dloglikeseq += Math.log(dscale); 

           if (bscaleemissions)
	   {
	      dloglikeseq += emissionproducts_scale[traindataObservedIndex_nseq[0]];
	   }

	   //forward step
	   //int numtime_nseq = numtime[nseq];
	   for (int nt = 1; nt < numtime_nseq; nt++)
	   {
	      //the actual observed combination at position t	        
	      double[] alpha_ntm1 = alpha[nt-1];
	      alpha_nt = alpha[nt];

	      
	      dscale = 0;
	      emissionproducts_nobserveindex = emissionproducts[traindataObservedIndex_nseq[nt]];
	      for (int ns = 0; ns < numstates; ns++)
	      {
	         //stores the emission product for each location on the chromosome		   

		 int transitionprobsnumCol_ns = transitionprobsnumCol[ns];
	         int[] transitionprobsindexCol_ns = transitionprobsindexCol[ns];
	         double[] coltransitionprobs_ns = coltransitionprobs[ns];

	         double dtempsum = 0;
                 if (transitionprobsnumCol_ns < nsparsecutoff)
		 {
	            //number of transitions is sparse enough worth going through the extra redirection
	            for (int nj = 0; nj < transitionprobsnumCol_ns; nj++)
	            {
	               //for each next state computing inner sum of all previous alpha and the transition probability
	               //for all non-zero transitions into the state
	               int nmappedindex = transitionprobsindexCol_ns[nj];
	               dtempsum += coltransitionprobs_ns[nmappedindex]*alpha_ntm1[nmappedindex];
	            }
		 }
		 else
	         {
	            //avoid the redirect and multiply by 0
                    for (int nj = 0; nj < numstates; nj++)
	            {
	               dtempsum += coltransitionprobs_ns[nj]*alpha_ntm1[nj];
	            }
		 }
	         //multiply the transition sum by the emission probability
	         double dalphaval = dtempsum*emissionproducts_nobserveindex[ns];
	         alpha_nt[ns] = dalphaval;
	         //System.out.println(ns+"\t"+alpha_nt[ns]+"\t"+dtempsum+"\t"+emissionproducts_nobserveindex[ns]);
	         dscale += dalphaval;
	      }

	      //rescaling alpha
	      scale[nt] = dscale;
              //scale_t(s)=P(o_0,...,o_t|lambda) summed over all states

	      if (bscalebeta)
	      {
	         if (dscale == 0)
	         {
		    for (int ns = 0; ns < numstates; ns++)
		    {
		       //added
		       if ((alpha_nt[ns] < EPSILONSTATE)&&(emissionproducts_nobserveindex[ns]>0))//(ns < numstates-1))) 
		       {
		          alpha_nt[ns] = EPSILONSTATE;
			  dscale += EPSILONSTATE;
		       }
		    }
		    scale[nt] = dscale;
		 }

	         for (int ns = 0; ns < numstates; ns++)
	         {
                    alpha_nt[ns] /= dscale;

  		    if ((alpha_nt[ns] < EPSILONSTATE)&&(emissionproducts_nobserveindex[ns]>0))//(ns < numstates-1))) 
		    {
		       alpha_nt[ns] = EPSILONSTATE;
		    }
		 }
	      }
	      else
	      {
	         for (int ns = 0; ns < numstates; ns++)
	         {
                    alpha_nt[ns] /= dscale;
		 }
	      }
      	      dloglikeseq += Math.log(dscale);

              if (bscaleemissions)
	      {
	         dloglikeseq += emissionproducts_scale[traindataObservedIndex_nseq[nt]];
	      }

	   }
	    
	   //backward step
	   //beta_t(s)=P(o_t+1,...,o_T|x_t=s,lambda)
           int nlastindex = numtime_nseq-1;
           double dinitval;

	   if (bscalebeta)
	   {
	       dinitval = 1.0/numstates;
	   }
	   else
	   {
               dinitval = 1.0/scale[nlastindex];
	   }

           for (int ns = 0; ns < numstates; ns++)
	   {
	      beta_ntp1[ns] = dinitval;
	   }

	   double ddenom = 0;	       

	   //gamma_nt - P(x=S| o_0,...,o_t)
	   //P(o_t+1,...,o_T|x_t=s,lambda) * P(o_0,...,o_t,xt=s|lambda)
	   alpha_nt = alpha[nlastindex];
	   for (int ns = 0; ns < gamma_nt.length; ns++)
           {
	       double dval = alpha_nt[ns]*beta_ntp1[ns];
	       ddenom += dval;
	       gamma_nt[ns] = dval;
	   }

	   if (ddenom > 0)
	   {
              for (int ns = 0; ns < gamma_nt.length; ns++)
	      {
	         gamma_nt[ns] /=ddenom;
	      }
	   }

	   double[] gammaObservedSum_combo_nt = gammaObservedSum[traindataObservedIndex_nseq[nlastindex]];
		
	   for (int ns = 0; ns < numstates; ns++)
	   { 
              //first sum gamma over all common signatures		     
	      //updates probability of observing the signature when in the state
	      gammaObservedSum_combo_nt[ns] += gamma_nt[ns];
	   }

	   for (int nt = nlastindex - 1; nt >= 0; nt--)
	   {
	      int ntp1 = (nt+1);
		   
	      double[] emissionproducts_combo_ntp1 = emissionproducts[traindataObservedIndex_nseq[ntp1]];

	      for (int ns = 0; ns < numstates; ns++)
	      {
	         tempproductbetaemiss[ns] = beta_ntp1[ns]*emissionproducts_combo_ntp1[ns];
	      }

	      if (bscaleemissions)
	      {
		  //adding here to help numerical stability
		  double dmaxval = 0;
		  for (int ns = 0; ns < numstates; ns++)
		  {
		      if (tempproductbetaemiss[ns] > dmaxval)
		      {
			  dmaxval = tempproductbetaemiss[ns];
		      }
		  }
                  for (int ns = 0; ns < numstates; ns++)
		  {
		      tempproductbetaemiss[ns] = tempproductbetaemiss[ns]/dmaxval;
		  }
	      }

	      double dsumbeta = 0;
	      //double dscaleinv = 1.0/scale[nt];
	      double dscale_nt = scale[nt];
              //scale_t(s)=P(o_0,...,o_t|lambda) summed over all states
	      for (int ni = 0; ni < numstates; ni++)
	      {
	         double dtempsum = 0;
	         int[] transitionprobsindex_ni =  transitionprobsindex[ni];
	         double[] transitionprobs_ni = transitionprobs[ni];
	         int transitionprobsnum_ni = transitionprobsnum[ni];

                 if (transitionprobsnum_ni < nsparsecutoff)
	         {
		    //sparse enought to pay the indirection penalty
		    for (int nj = 0; nj < transitionprobsnum_ni; nj++)
		    {
		       //for each state summing over transition probability to state j, emission probablity in j at next step
		       //and probability of observing the remaining sequence
		       int nmappedindexouter = transitionprobsindex_ni[nj];
		       dtempsum += transitionprobs_ni[nmappedindexouter]*tempproductbetaemiss[nmappedindexouter];
		    }
		 }
		 else
		 {
		    //not trying to exploit sparsity here
                    for (int nj = 0; nj < numstates; nj++)
		    {
		       //for each state summing over transition probability to state j, emission probablity in j at next step
		       //and probability of observing the remaining sequence
		       dtempsum += transitionprobs_ni[nj]*tempproductbetaemiss[nj];
		    }
		 }

		 if (bscalebeta)
		 {
	            beta_nt[ni] = dtempsum;
		    dsumbeta += dtempsum;
		 }
	         else
		 {
		    double dratio = dtempsum/dscale_nt;

		    if (dratio > Double.MAX_VALUE)
		    {
		       beta_nt[ni] = Double.MAX_VALUE;// dtempsum/dscale_nt;
		    }
		    else
		    {
		       beta_nt[ni] = dratio;// dtempsum/dscale_nt;
		    }
		 }
	      }

	      if (bscalebeta)
	      {
                 for (int ni = 0; ni < numstates; ni++)
	         {
	            beta_nt[ni]/= dsumbeta;

		    if (beta_nt[ni] < EPSILONSTATE)//&&(!bdummy))// || (ni < numstates-1))) 
		    {
		       beta_nt[ni] = EPSILONSTATE;
		    }
		 }
	      }		
		
	      ddenom = 0;
		
	      alpha_nt = alpha[nt];	    

	      //gamma_nt - P(x=S| o_0,...,o_t)
	      //P(o_t+1,...,o_T|x_t=s,lambda) * P(o_0,...,o_t,xt=s|lambda)

	      for (int ns = 0; ns < gamma_nt.length; ns++)
	      {
	         double dval = alpha_nt[ns]*beta_nt[ns];
	         ddenom += dval;
	         gamma_nt[ns] = dval;
	      }

	      if (ddenom > 0)
	      {
	         for (int ns = 0; ns < gamma_nt.length; ns++)
	         {
	            gamma_nt[ns] /= ddenom;
		 }
	      }

              gammaObservedSum_combo_nt = gammaObservedSum[traindataObservedIndex_nseq[nt]];

              for (int ns = 0; ns < numstates; ns++)
	      {
	         //first sum gamma over all common signatures
		 //updates probability of observing the signature when in the state
	         gammaObservedSum_combo_nt[ns] += gamma_nt[ns];
	      }

	      double dsum = 0;
		  
              //this compues the numerator portion
	      for (int ni = 0; ni < numstates; ni++)
	      {
	         double[] sumforsxi_ni = sumforsxi[ni]; //computing expected number of transition from state i
	         int[] transitionprobsindex_ni = transitionprobsindex[ni]; //indicies of non-zero transitions from state i
	         double[] transitionprobs_ni = transitionprobs[ni]; //probability of transitions from state i
	         int ntransitionprobsnum_ni = transitionprobsnum[ni]; //number of non-zero transitions from state i
	         double dalpha_nt_ni = alpha_nt[ni]; 
	         //sxi is P(q_t = S_i, q_(t+1) = S_j | O)

                 if (ntransitionprobsnum_ni < nsparsecutofflooser)
      	         {
		    //enough 0 transitionto use sparsity here
		    //looser cut off since the indirection is less of the total time
	            for (int nj = 0; nj < ntransitionprobsnum_ni; nj++)
		    {
		       int nmappedindex = transitionprobsindex_ni[nj];
		       //computes transition probability from state i to j
		       double dtempval = transitionprobs_ni[nmappedindex] *dalpha_nt_ni*tempproductbetaemiss[nmappedindex];
                       dsum += dtempval;
	      	       sumforsxi_ni[nmappedindex] = dtempval;
		    }
		 }
	         else
	         {
		    for (int nj = 0; nj < numstates; nj++)
		    {
		       //computes transition probability from state i to j
		       double dtempval = transitionprobs_ni[nj]*dalpha_nt_ni*tempproductbetaemiss[nj];
		       dsum += dtempval;
		       sumforsxi_ni[nj] = dtempval;
		    }
		 }
	      }    	 
		   
	      if (dsum > 0)
	      {
	         //normalizing the numerator by the sum of the denominator and updating this iterations value for it
	         for (int ni = 0; ni < numstates; ni++)
	         {
		    int[] transitionprobsindex_ni = transitionprobsindex[ni];
		    double[] sumforsxi_ni = sumforsxi[ni];
		    double[] sxi_nseq_ni = sxi_nseq[ni];
		    int ntransitionprobsnum_ni = transitionprobsnum[ni];

                    if (ntransitionprobsnum_ni < nsparsecutoff)
	            {
		       //guessing sparse enough to avoid the indirections
		       for (int nj = 0; nj < ntransitionprobsnum_ni; nj++)
		       {
	      	          int nmappedindex = transitionprobsindex_ni[nj];
		          sxi_nseq_ni[nmappedindex] += sumforsxi_ni[nmappedindex]/dsum;
		       }
		    }
	            else
	            {
                       for (int nj = 0; nj < numstates; nj++)
	      	       {
		          sxi_nseq_ni[nj] += sumforsxi_ni[nj]/dsum;
		       }
		    }
		 }
	      }   	 
	      beta_ntp1 = beta_nt;  //updating beta_ntp1 
	   }

	   //double[] gammainitstore_nseq = gammainitstore[nseq];
	   for (int ns = 0; ns < numstates; ns++)
	   {
	      //storing the initial gamma from this iteration
	      gammainitstore_nseq[ns] = gamma_nt[ns];
	   }

	  for (int nindex = 0; nindex < gammaObservedSum.length; nindex++)
	  {
	     //going through all the gamma sufficient statistic
	     if (traindataObservedSeqFlags_nseq[nindex])
	     {
	        //only update for those combinations that were observed on this sequnce
		//gets the observed combination and missing combination signatures
	        boolean[] traindataObservedValues_nindex = traindataObservedValues[nindex];
	        boolean[] traindataNotMissing_nindex = traindataNotMissing[nindex];
	           
		double[] gammaObservedSum_nindex = gammaObservedSum[nindex];

		for (int ns = 0; ns < numstates; ns++)
	        {
	           //going through each state
		   double[][] gammaksum_nseq_ns = gammaksum_nseq[ns];
	           double gammaObservedSum_nindex_ns = gammaObservedSum_nindex[ns];
		   for (int nmark = 0; nmark < numdatasets; nmark++)
	           {
		      //going through each mark
		      if (traindataNotMissing_nindex[nmark])
		      {
		         //only update non-missing
		         if (traindataObservedValues_nindex[nmark])
		         {
		            //updates the gamma sum for each mark when in state and observed 1
		            gammaksum_nseq_ns[nmark][1] += gammaObservedSum_nindex_ns;
			 }
		         else
		         {
		    	    //updates the gamma sum for each mark when in state and observed 0
		            gammaksum_nseq_ns[nmark][0] += gammaObservedSum_nindex_ns;
			 }
		      }
		   }
		}
	     }
	  }
          dloglikeA[nseq] = dloglikeseq;	
	}
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////




    /////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create a new thread.
    class NewThreadWithLoad implements Runnable 
    {
       
	//int[] traindataObservedIndex_nseq;
	//boolean[] traindataObservedSeqFlags_nseq;
	double[][][] gammaksum_nseq;
	double[][] sxi_nseq;
	int numtime_nseq;
	String chromfiles_nseq;
	double[] gammainitstore_nseq;
	double[][][] gammaObservedSum_Pool;
	double[][][] alpha_Pool;
	//double[][] emissionproducts;
	double[][] gamma_nt_Pool;
	double[][] coltransitionprobs;
	double[][] scale_Pool;
	double[][] beta_nt_Pool;
	double[][] beta_ntp1_Pool;
	double[][] tempproductbetaemiss_Pool;
	int nsparsecutoff;
	int nsparsecutofflooser;
	double[][][] sumforsxi_Pool;
	double[] dloglikeA;
	int nseq;

	int[][] traindataObservedIndex_Pool;
        boolean[][][] traindataObservedValues_Pool;
        boolean[][][] traindataNotMissing_Pool;
	double[][][] emissionproducts_Pool;
        double[][] emissionproducts_scale_Pool;


	NewThreadWithLoad(
                          String chromfiles_nseq,
                          int[][] traindataObservedIndex_Pool,
			  boolean[][][] traindataObservedValues_Pool,
			  boolean[][][] traindataNotMissing_Pool,
			  //boolean[] traindataObservedSeqFlags_nseq,
		  double[][][] gammaksum_nseq,
		  double[][] sxi_nseq,
		  int numtime_nseq,
		  double[] gammainitstore_nseq,
		  double[][][] gammaObservedSum_Pool,
		  double[][][] alpha_Pool,
		  double[][][] emissionproducts_Pool,
		  double[][] gamma_nt_Pool,
		  double[][] coltransitionprobs,
		  double[][] scale_Pool,
		  double[][] beta_nt_Pool,
		  double[][] beta_ntp1_Pool,
		  double[][] tempproductbetaemiss_Pool,
		  double[][][] sumforsxi_Pool,
		  int nsparsecutoff,
		  int nsparsecutofflooser,
		  double[] dloglikeA,
		  int nseq,
                  double[][] emissionproducts_scale_Pool)
        {
	    this.emissionproducts_Pool = emissionproducts_Pool;
            this.emissionproducts_scale_Pool = emissionproducts_scale_Pool;
	    this.chromfiles_nseq = chromfiles_nseq;
	    this.traindataObservedIndex_Pool = traindataObservedIndex_Pool;
	    this.traindataObservedValues_Pool = traindataObservedValues_Pool;
	    this.traindataNotMissing_Pool = traindataNotMissing_Pool;
	    //this.traindataObservedIndex_nseq = traindataObservedIndex_nseq;
	    //this.traindataObservedSeqFlags_nseq = traindataObservedSeqFlags_nseq;
	    this.gammaksum_nseq = gammaksum_nseq;
	    this.sxi_nseq = sxi_nseq;
	    this.numtime_nseq = numtime_nseq;
	    this.gammainitstore_nseq = gammainitstore_nseq;
	    this.gammaObservedSum_Pool = gammaObservedSum_Pool;
	    this.alpha_Pool = alpha_Pool;
	    //this.emissionproducts = emissionproducts;
	    this.gamma_nt_Pool = gamma_nt_Pool;
	    this.coltransitionprobs = coltransitionprobs;
	    this.scale_Pool = scale_Pool;
	    this.beta_nt_Pool = beta_nt_Pool;
	    this.beta_ntp1_Pool = beta_ntp1_Pool;
	    this.tempproductbetaemiss_Pool = tempproductbetaemiss_Pool;
	    this.nsparsecutoff = nsparsecutoff;
	    this.nsparsecutofflooser = nsparsecutofflooser;
	    this.sumforsxi_Pool = sumforsxi_Pool;
	    this.dloglikeA = dloglikeA;
	    this.nseq = nseq;
	}

	private int slotavailable()
	{
	   for (int ni = 0; ni < threadslots.length; ni++)
	   {
	      if (!threadslots[ni])
	      {
	         threadslots[ni] = true;
	         return ni;
	      }
	   }
	   return -1;	    
	}

	// This is the entry point for the second thread.
	public void run() 
        {
	    int nprocess;

	    synchronized (objlock)
	    {
	       while ((nprocess=slotavailable())==-1)
	       {
                  try
                  {
	             objlock.wait();
	          }
	          catch (InterruptedException ex) {}
	       }
	    }

            try
            {

  	       estep(
		  emissionproducts_Pool[nprocess],
		  traindataObservedIndex_Pool[nprocess],
		  traindataObservedValues_Pool[nprocess],
		  traindataNotMissing_Pool[nprocess],
                  gammaObservedSum_Pool[nprocess],
		  alpha_Pool[nprocess],
		  gamma_nt_Pool[nprocess],
		  scale_Pool[nprocess],
		  beta_nt_Pool[nprocess],
		  beta_ntp1_Pool[nprocess],
		  tempproductbetaemiss_Pool[nprocess],
		  sumforsxi_Pool[nprocess],
                  emissionproducts_scale_Pool[nprocess]);
	    }
	    catch (IOException ioex)
	    {
		ioex.printStackTrace(System.out);
	    }
	    
	    synchronized(objlock)
	    {
		nlaunched--;
		threadslots[nprocess] = false;
		objlock.notifyAll();
	    }
	}

        public void estep(
		  double[][] emissionproducts,
		  int[] traindataObservedIndex,
		  boolean[][] traindataObservedValues,
		  boolean[][] traindataNotMissing,
		  double[][] gammaObservedSum,
		  double[][] alpha,
		  double[] gamma_nt,
		  double[] scale,
		  double[] beta_nt,
		  double[] beta_ntp1,
		  double[] tempproductbetaemiss,
		  double[][] sumforsxi,
                  double[] emissionproducts_scale) throws IOException
	{ 
	      //going through each sequence

	      ///////////////////////////////////////////////////////////
	      //loading the data just for the current chromosome 

	      HashMap hmObserved = new HashMap(); //maps an observation string to an index and set of flags
 
    	      int nobserved = 0;

	      if (ChromHMM.BVERBOSE)
	      {
                 System.out.println("reading\t"+szinputdir+" "+chromfiles_nseq);
	      }

	      BufferedReader br = Util.getBufferedReader(szinputdir+"/"+chromfiles_nseq);
              String szLine; 
              br.readLine(); //first line tells cell type and chromosome
	      br.readLine();//flush mark header
	      ArrayList aldata = new ArrayList();
	      while ((szLine = br.readLine())!=null)
	      {
	         StringTokenizer st = new StringTokenizer(szLine,"\t ");
		 StringBuffer sb = new StringBuffer();
		
		 for (int ncol = 0; ncol < numdatasets; ncol++)
		 {
		    if (!st.hasMoreTokens())
		    {
	      	       throw new IllegalArgumentException("Found line without "+numdatasets+" values in file "+chromfiles_nseq);
		    }

		    String sztoken = st.nextToken();//.trim();
		    
		    if (sztoken.equals("0"))
		    {
		       sb.append("0");
		    }
		    else if (sztoken.equals("1"))
		    {
		       sb.append("1");
		    }
		    else if (sztoken.equals("2"))
		    {
		       //this means missing
	               sb.append("2");
		    }
		    else
	            {
	               throw new IllegalArgumentException("Unrecognized value "+sztoken+" found in "+szinputdir+"/"+chromfiles_nseq);
	            }
		 }
		 aldata.add(sb.toString());
	      }
	      br.close();
	      
	      int nsize = aldata.size();
	      //traindataObservedIndex[nfile] = new int[nsize];
	      //int[] traindataObservedIndex_nfile = traindataObservedIndex[nfile];

	      for (int nrow = 0; nrow < nsize; nrow++)
	      {
	         BigInteger theBigInteger = new BigInteger((String) aldata.get(nrow),3);
		 Integer theObservedInt  = (Integer) hmObserved.get(theBigInteger);
	         //boolean[] flagA;

	         if (theObservedInt == null)
       	         {
		    //System.out.println(szmappingbyte.length());
		    //storing a mapping from observed byte string to an integer index in alFlags and alObserved
		     hmObserved.put(theBigInteger, Integer.valueOf(nobserved));

		    //saving this observed index
		    traindataObservedIndex[nrow] = nobserved;

		    //increments the number of observed combinations of marks
		    nobserved++;
		 }
		 else
		 {
		    //storing the index of the flags associated with this row 
		     traindataObservedIndex[nrow] = ((Integer) theObservedInt).intValue();
		 }
	      }	  	 
	   


	     Iterator hmObservedIterator = hmObserved.entrySet().iterator();
	     while (hmObservedIterator.hasNext())
	     {
	        Map.Entry pairs = (Map.Entry) hmObservedIterator.next();
	        BigInteger theBigInteger = (BigInteger) pairs.getKey();
	        String szmapping = theBigInteger.toString(3);  //getting back the mapping string

	        //ObservedRec theObservedRec = (ObservedRec) pairs.getValue();

	        int ncurrindex = ((Integer) pairs.getValue()).intValue();// theObservedRec.nobserved;//this is an index on which obervation combination it is

	        boolean[] traindataObservedValues_ncurrindex = traindataObservedValues[ncurrindex];
	        boolean[] traindataNotMissing_ncurrindex = traindataNotMissing[ncurrindex]; 
	    
	        //if the mapping string is less than the number of data sets then 
	        //there are leading 0's will set for leading 0's not missing and absent
	        int numch = szmapping.length();
	        int numleading0 = numdatasets - numch;
	        for (int nj = 0; nj < numleading0; nj++)
	        {
		   traindataObservedValues_ncurrindex[nj] = false;
		   traindataNotMissing_ncurrindex[nj] = true;
		}

	        int nmappedindex = numleading0; //starting from the leading 0 position
	        for (int nj = 0; nj < numch; nj++)
	        {
	           char ch = szmapping.charAt(nj);

	           if (ch == '0')
	           {
		      traindataObservedValues_ncurrindex[nmappedindex] = false;
		      traindataNotMissing_ncurrindex[nmappedindex] = true;
		   }
	           else if (ch=='1')
	           {
		      traindataObservedValues_ncurrindex[nmappedindex] = true;
		      traindataNotMissing_ncurrindex[nmappedindex] = true;
		   }
	           else
	           {
		      //missing data
		      traindataObservedValues_ncurrindex[nmappedindex] = false;
		      traindataNotMissing_ncurrindex[nmappedindex] = false;
		   }
	           nmappedindex++;
		}
	     }     

	     if (bscaleemissions)
	     {
	        for (int ni = 0; ni < nobserved; ni++)
		{
		   //going through each combination of marks
		   //if (traindataObservedSeqFlags_nseq[ni])
		   //{
		   //this signature of marks is observed on the current chromosome so
		   //updating its emission probabilities
		   double[] emissionproducts_ni = emissionproducts[ni];
		   boolean[] traindataObservedValues_ni = traindataObservedValues[ni];
		   boolean[] traindataNotMissing_ni = traindataNotMissing[ni];

		   for (int ns = 0; ns < numstates; ns++)
		   {
		      emissionproducts_ni[ns] = 1;
		   }
		   emissionproducts_scale[ni] = 0;

		   for (int nmod = 0; nmod < numdatasets; nmod++)
		   {
		      for (int ns = 0; ns < numstates; ns++)
		      {
		         if (traindataNotMissing_ni[nmod])
			 {
			    //we are include this marks emission probability
			    if (traindataObservedValues_ni[nmod])
			    {
				emissionproducts_ni[ns] *= emissionprobs[ns][nmod][1];
			    }
			    else
			    {
				emissionproducts_ni[ns] *= emissionprobs[ns][nmod][0];
			    }
			 }
					 // otherwise treated as missing omitting from product
		      }

		      
		      double dmaxval = 0;
		      for (int ns = 0; ns < numstates; ns++)
		      {
		         if (emissionproducts_ni[ns] > dmaxval)
			 {
			    dmaxval = emissionproducts_ni[ns];
			 }
		      }

		      if (dmaxval <=0)///EPSILONEMISSIONS)
		      {
		         for (int ns = 0; ns < numstates; ns++)
		         {
		            emissionproducts_ni[ns] = 1;
		         }
			 //approximating log-likelihood
	      	         emissionproducts_scale[ni] += Math.log(EPSILONEMISSIONS);
		      }
		      else
		      {
		         for (int ns = 0; ns < numstates; ns++)
			 {
		            emissionproducts_ni[ns]/= dmaxval;
			 }
			 emissionproducts_scale[ni] += Math.log(dmaxval);			  
		      }		      	      
		   }		
		}
	     }
	     else
	     {
	        for (int ni = 0; ni < nobserved; ni++)
	        {
	           //going through each combination of marks
	           //if (traindataObservedSeqFlags_nseq[ni])
		   //{
	           //this signature of marks is observed on the current chromosome so
       	           //updating its emission probabilities
       	           double[] emissionproducts_ni = emissionproducts[ni];
      	           boolean[] traindataObservedValues_ni = traindataObservedValues[ni];
       	           boolean[] traindataNotMissing_ni = traindataNotMissing[ni];		  

	           boolean ballzero = true;

	           for (int ns = 0; ns < numstates; ns++)
	           {
	              double dproduct = 1;
	              double[][] emissionprobs_ns = emissionprobs[ns];

		      for (int nmod = 0; nmod < numdatasets; nmod++)
	              {
		         if (traindataNotMissing_ni[nmod])
	      	         {
		            //we are include this marks emission probability
		            if (traindataObservedValues_ni[nmod])
		            {
		      	       //System.out.println("positive\t"+ns+"\t"+nmod+"\t1\t"+emissionprobs_ns[nmod][1]);
		               dproduct *= emissionprobs_ns[nmod][1];
		            }
		            else 
		            {
		      	       ///System.out.println("negative\t"+ns+"\t"+nmod+"\t0\t"+emissionprobs_ns[nmod][0]);
		               dproduct *= emissionprobs_ns[nmod][0];
			    }
			 }
			   // otherwise treated as missing omitting from product
		      }
		      //System.out.println(ns+"\t"+dproduct);
		      emissionproducts_ni[ns] = dproduct;

		      if (dproduct >= EPSILONEMISSIONS)
		      {
	                 ballzero = false;
		      }
		   }

		   if (ballzero)
	           {
	              for (int ns = 0; ns < numstates; ns++)
		      {
		         emissionproducts_ni[ns] = EPSILONEMISSIONS;
		      }
		   }
		}
	     }



	   double dloglikeseq = 0;

	   for (int ns = 0; ns < gammaksum_nseq.length; ns++)
	   {
	      //resetting the gamma sufficient statistics in the current sequence
	      double[][] gammaksum_nseq_ns = gammaksum_nseq[ns];
	      for (int nmark = 0; nmark < gammaksum_nseq_ns.length; nmark++)
	      {
	         for (int nbucket = 0; nbucket < numbuckets; nbucket++)
	         {
	            gammaksum_nseq_ns[nmark][nbucket] = 0;
	         }
	      }
           }

	   for (int ni = 0; ni < sxi_nseq.length; ni++)
	   {
	      //reseeting the sxi sufficient statistics in the current sequence
	      double[] sxi_nseq_ni = sxi_nseq[ni];
	      for (int nj = 0; nj < sxi_nseq_ni.length; nj++)
	      {
                 sxi_nseq_ni[nj] = 0;
	      }
	   }

	   //gammaObservedSum stores the weight for each combination of marks in each stat
	   //for (int ncombo = 0; ncombo < gammaObservedSum.length; ncombo++)
	   for (int ncombo = 0; ncombo < nobserved; ncombo++)
	   {
	      //resetting that to 0
	      double[] gammaObservedSum_ncombo = gammaObservedSum[ncombo];
	      for (int ns = 0; ns < gammaObservedSum_ncombo.length; ns++)
	      {
	         gammaObservedSum_ncombo[ns] = 0;
	      }
	   }

	   //initial probability in state s is initial probability times emission probability at first position
	   double[] alpha_nt = alpha[0];
           double dscale = 0;
	   double[] emissionproducts_nobserveindex =emissionproducts[traindataObservedIndex[0]];
 	   for (int ns = 0; ns < numstates; ns++)
           {
	      alpha_nt[ns] = probinit[ns] * emissionproducts_nobserveindex[ns];
	      //System.out.println(ns+"\t"+alpha_nt[ns]);
	      dscale += alpha_nt[ns];
	   }
	   scale[0] = dscale;

	   //alpha_t(s)=P(o_0,...,o_t,x_t=s|lambda)
	   //converts the alpha terms to probabilities
	   if (bscalebeta)
           {
	      if (dscale == 0)
	      {
	         for (int ns = 0; ns < numstates; ns++)
		 {
		    //added
		    if ((alpha_nt[ns] < EPSILONSTATE)&&(emissionproducts_nobserveindex[ns]>0))//(ns < numstates-1))) 
		    {
		       alpha_nt[ns] = EPSILONSTATE;
	       	       dscale += EPSILONSTATE;
	            }
	         }
	         scale[0] = dscale;
	      }

	      for (int ns = 0; ns < numstates; ns++)
	      {
                 alpha_nt[ns] /= dscale;

  		 if ((alpha_nt[ns] < EPSILONSTATE)&&(emissionproducts_nobserveindex[ns]>0))//(ns < numstates-1))) 
	         {
	            alpha_nt[ns] = EPSILONSTATE;
		 }
	      }
	   }
	   else
	   {
	      for (int ns = 0; ns < numstates; ns++)
	      {
                 alpha_nt[ns] /= dscale;
	      }
	   }
	     


	   dloglikeseq += Math.log(dscale); 

           if (bscaleemissions)
	   {
	      dloglikeseq += emissionproducts_scale[traindataObservedIndex[0]];
	   }

	   //forward step
	   //int numtime_nseq = numtime[nseq];
	   for (int nt = 1; nt < numtime_nseq; nt++)
	   {
	      //the actual observed combination at position t	        
	      double[] alpha_ntm1 = alpha[nt-1];
	      alpha_nt = alpha[nt];
	      
	      dscale = 0;
	      emissionproducts_nobserveindex = emissionproducts[traindataObservedIndex[nt]];

	      for (int ns = 0; ns < numstates; ns++)
	      {
	         //stores the emission product for each location on the chromosome		   

		 int transitionprobsnumCol_ns = transitionprobsnumCol[ns];
	         int[] transitionprobsindexCol_ns = transitionprobsindexCol[ns];
	         double[] coltransitionprobs_ns = coltransitionprobs[ns];

	         double dtempsum = 0;
                 if (transitionprobsnumCol_ns < nsparsecutoff)
		 {
	            //number of transitions is sparse enough worth going through the extra redirection
	            for (int nj = 0; nj < transitionprobsnumCol_ns; nj++)
	            {
	               //for each next state computing inner sum of all previous alpha and the transition probability
	               //for all non-zero transitions into the state
	               int nmappedindex = transitionprobsindexCol_ns[nj];
	               dtempsum += coltransitionprobs_ns[nmappedindex]*alpha_ntm1[nmappedindex];
	            }
		 }
		 else
	         {
	            //avoid the redirect and multiply by 0
                    for (int nj = 0; nj < numstates; nj++)
	            {
	               dtempsum += coltransitionprobs_ns[nj]*alpha_ntm1[nj];
	            }
		 }
	         //multiply the transition sum by the emission probability
	         double dalphaval = dtempsum*emissionproducts_nobserveindex[ns];
	         alpha_nt[ns] = dalphaval;
	         //System.out.println(ns+"\t"+alpha_nt[ns]+"\t"+dtempsum+"\t"+emissionproducts_nobserveindex[ns]);
	         dscale += dalphaval;
	      }


	      //rescaling alpha
	      scale[nt] = dscale;
              //scale_t(s)=P(o_0,...,o_t|lambda) summed over all states

	      if (bscalebeta)
	      {
	         if (dscale == 0)
	         {
		    for (int ns = 0; ns < numstates; ns++)
		    {
		       //added
		       if ((alpha_nt[ns] < EPSILONSTATE)&&(emissionproducts_nobserveindex[ns]>0))//(ns < numstates-1))) 
		       {
		          alpha_nt[ns] = EPSILONSTATE;
			  dscale += EPSILONSTATE;
		       }
		    }
		    scale[nt] = dscale;
		 }

		 //		 else
		 //{
	         for (int ns = 0; ns < numstates; ns++)
	         {
                    alpha_nt[ns] /= dscale;

		    if ((alpha_nt[ns] < EPSILONSTATE)&&(emissionproducts_nobserveindex[ns]>0))//(ns < numstates-1))) 
		    {
	               alpha_nt[ns] = EPSILONSTATE;
		    }
		 }		 
		 //}
	      }
      	      else
	      {
	         for (int ns = 0; ns < numstates; ns++)
	         {
                    alpha_nt[ns] /= dscale;
		 }
	      }
	      //for (int ns = 0; ns < numstates; ns++)
	      //{
              //   alpha_nt[ns] /= dscale;
              //}

      	      dloglikeseq += Math.log(dscale);

              if (bscaleemissions)
	      {
	         dloglikeseq += emissionproducts_scale[traindataObservedIndex[nt]];
	      }
	   }
	    
	   //backward step
	   //beta_t(s)=P(o_t+1,...,o_T|x_t=s,lambda)
           int nlastindex = numtime_nseq-1;
	   double dinitval;
	   if (bscalebeta)
	   {
	      dinitval = 1.0/numstates;
	   }
	   else
	   {
              dinitval = 1.0/scale[nlastindex];
	   }

           for (int ns = 0; ns < numstates; ns++)
	   {
	      beta_ntp1[ns] = dinitval;
	   }
	   //double dinitval = 1.0/scale[nlastindex];
	   //for (int ns = 0; ns < numstates; ns++)
	   //{
	   // beta_ntp1[ns] = dinitval;
	   //}

	   double ddenom = 0;	       

	   //gamma_nt - P(x=S| o_0,...,o_t)
	   //P(o_t+1,...,o_T|x_t=s,lambda) * P(o_0,...,o_t,xt=s|lambda)
	   alpha_nt = alpha[nlastindex];
	   for (int ns = 0; ns < gamma_nt.length; ns++)
           {
	       double dval = alpha_nt[ns]*beta_ntp1[ns];
	       ddenom += dval;
	       gamma_nt[ns] = dval;
	   }

	   if (ddenom > 0)
	   {
              for (int ns = 0; ns < gamma_nt.length; ns++)
	      {
	         gamma_nt[ns] /=ddenom;
	      }
	   }

	   double[] gammaObservedSum_combo_nt = gammaObservedSum[traindataObservedIndex[nlastindex]];
		
	   for (int ns = 0; ns < numstates; ns++)
	   { 
              //first sum gamma over all common signatures		     
	      //updates probability of observing the signature when in the state
	      gammaObservedSum_combo_nt[ns] += gamma_nt[ns];
	   }

	   for (int nt = nlastindex - 1; nt >= 0; nt--)
	   {
	      int ntp1 = (nt+1);
		   
	      double[] emissionproducts_combo_ntp1 = emissionproducts[traindataObservedIndex[ntp1]];

	      for (int ns = 0; ns < numstates; ns++)
	      {
	         tempproductbetaemiss[ns] = beta_ntp1[ns]*emissionproducts_combo_ntp1[ns];
	      }

	      if (bscaleemissions)
	      {
		  //adding here to help numerical stability
		  double dmaxval = 0;
		  for (int ns = 0; ns < numstates; ns++)
		  {
		      if (tempproductbetaemiss[ns] > dmaxval)
		      {
			  dmaxval = tempproductbetaemiss[ns];
		      }
		  }
                  for (int ns = 0; ns < numstates; ns++)
		  {
		      tempproductbetaemiss[ns] = tempproductbetaemiss[ns]/dmaxval;
		  }
	      }

	      double dsumbeta = 0;

	      //double dscaleinv = 1.0/scale[nt];
	      double dscale_nt = scale[nt];
              //scale_t(s)=P(o_0,...,o_t|lambda) summed over all states
	      for (int ni = 0; ni < numstates; ni++)
	      {
	         double dtempsum = 0;
	         int[] transitionprobsindex_ni =  transitionprobsindex[ni];
	         double[] transitionprobs_ni = transitionprobs[ni];
	         int transitionprobsnum_ni = transitionprobsnum[ni];

                 if (transitionprobsnum_ni < nsparsecutoff)
	         {
		    //sparse enought to pay the indirection penalty
		    for (int nj = 0; nj < transitionprobsnum_ni; nj++)
		    {
		       //for each state summing over transition probability to state j, emission probablity in j at next step
		       //and probability of observing the remaining sequence
		       int nmappedindexouter = transitionprobsindex_ni[nj];
		       dtempsum += transitionprobs_ni[nmappedindexouter]*tempproductbetaemiss[nmappedindexouter];
		    }
		 }
		 else
		 {
		    //not trying to exploit sparsity here
                    for (int nj = 0; nj < numstates; nj++)
		    {
		       //for each state summing over transition probability to state j, emission probablity in j at next step
		       //and probability of observing the remaining sequence
		       dtempsum += transitionprobs_ni[nj]*tempproductbetaemiss[nj];
		    }
		 }

		 if (bscalebeta)
		 {
		    beta_nt[ni] = dtempsum;
		    dsumbeta += dtempsum;
		 }
	         else
		 {
	            double dratio = dtempsum/dscale_nt;
		    if (dratio > Double.MAX_VALUE)
		    {
		       beta_nt[ni] = Double.MAX_VALUE;//dtempsum/dscale_nt;
		    }
		    else
		    {
		       beta_nt[ni] = dratio;
		    }
		 }
		 //double dratio = dtempsum/dscale_nt;

		 //if (dratio > Double.MAX_VALUE)
		 //{
		 //    beta_nt[ni] = Double.MAX_VALUE;// dtempsum/dscale_nt;
		 //}
		 //else
		 //{
		 //    beta_nt[ni] = dratio;//dtempsum/dscale_nt;
		 //}
	      }

	      if (bscalebeta)
	      {
                 for (int ni = 0; ni < numstates; ni++)
	         {
	            beta_nt[ni]/= dsumbeta;

		    if (beta_nt[ni] < EPSILONSTATE)//&&(!bdummy))// || (ni < numstates-1))) 
		    {
		       beta_nt[ni] = EPSILONSTATE;
		    }
		 }
	      }		
		
	      ddenom = 0;
		
	      alpha_nt = alpha[nt];	    

	      //gamma_nt - P(x=S| o_0,...,o_t)
	      //P(o_t+1,...,o_T|x_t=s,lambda) * P(o_0,...,o_t,xt=s|lambda)

	      for (int ns = 0; ns < gamma_nt.length; ns++)
	      {
	         double dval = alpha_nt[ns]*beta_nt[ns];
	         ddenom += dval;
	         gamma_nt[ns] = dval;
	      }

	      if (ddenom > 0)
	      {
	         for (int ns = 0; ns < gamma_nt.length; ns++)
	         {
	            gamma_nt[ns] /= ddenom;
		 }
	      }

              gammaObservedSum_combo_nt = gammaObservedSum[traindataObservedIndex[nt]];

              for (int ns = 0; ns < numstates; ns++)
	      {
	         //first sum gamma over all common signatures
		 //updates probability of observing the signature when in the state
	         gammaObservedSum_combo_nt[ns] += gamma_nt[ns];
	      }

	      double dsum = 0;
		  
              //this compues the numerator portion
	      for (int ni = 0; ni < numstates; ni++)
	      {
	         double[] sumforsxi_ni = sumforsxi[ni]; //computing expected number of transition from state i
	         int[] transitionprobsindex_ni = transitionprobsindex[ni]; //indicies of non-zero transitions from state i
	         double[] transitionprobs_ni = transitionprobs[ni]; //probability of transitions from state i
	         int ntransitionprobsnum_ni = transitionprobsnum[ni]; //number of non-zero transitions from state i
	         double dalpha_nt_ni = alpha_nt[ni]; 
	         //sxi is P(q_t = S_i, q_(t+1) = S_j | O)

                 if (ntransitionprobsnum_ni < nsparsecutofflooser)
      	         {
		    //enough 0 transitionto use sparsity here
		    //looser cut off since the indirection is less of the total time
	            for (int nj = 0; nj < ntransitionprobsnum_ni; nj++)
		    {
		       int nmappedindex = transitionprobsindex_ni[nj];
		       //computes transition probability from state i to j
		       double dtempval = transitionprobs_ni[nmappedindex] *dalpha_nt_ni*tempproductbetaemiss[nmappedindex];
                       dsum += dtempval;

	      	       sumforsxi_ni[nmappedindex] = dtempval;
		    }
		 }
	         else
	         {
		    for (int nj = 0; nj < numstates; nj++)
		    {
		       //computes transition probability from state i to j
		       double dtempval = transitionprobs_ni[nj]*dalpha_nt_ni*tempproductbetaemiss[nj];
		       dsum += dtempval;
		       sumforsxi_ni[nj] = dtempval;
		    }
		 }
	      }

	      
	      if (dsum > 0)
	      {
		 //added to avoid numerical errors
	         //normalizing the numerator by the sum of the denominator and updating this iterations value for it
	         for (int ni = 0; ni < numstates; ni++)
	         {
		    int[] transitionprobsindex_ni = transitionprobsindex[ni];
		    double[] sumforsxi_ni = sumforsxi[ni];
		    double[] sxi_nseq_ni = sxi_nseq[ni];
		    int ntransitionprobsnum_ni = transitionprobsnum[ni];

                    if (ntransitionprobsnum_ni < nsparsecutoff)
	            {
		       //guessing sparse enough to avoid the indirections
		       for (int nj = 0; nj < ntransitionprobsnum_ni; nj++)
		       {
	      	          int nmappedindex = transitionprobsindex_ni[nj];
		          sxi_nseq_ni[nmappedindex] += sumforsxi_ni[nmappedindex]/dsum;
		       }
		    }
	            else
	            {
                       for (int nj = 0; nj < numstates; nj++)
	      	       {
		          sxi_nseq_ni[nj] += sumforsxi_ni[nj]/dsum;
		       }
		    }
		 }
	      }	 
	      beta_ntp1 = beta_nt;  //updating beta_ntp1 
	   }

	   //double[] gammainitstore_nseq = gammainitstore[nseq];
	   for (int ns = 0; ns < numstates; ns++)
	   {
	      //storing the initial gamma from this iteration
	      gammainitstore_nseq[ns] = gamma_nt[ns];
	   }

	   //for (int nindex = 0; nindex < gammaObservedSum.length; nindex++)
	  for (int nindex = 0; nindex < nobserved; nindex++)
	  {
	     //going through all the gamma sufficient statistic
	     //if (traindataObservedSeqFlags_nseq[nindex])
	     {
	        //only update for those combinations that were observed on this sequence
		//gets the observed combination and missing combination signatures
	        boolean[] traindataObservedValues_nindex = traindataObservedValues[nindex];
	        boolean[] traindataNotMissing_nindex = traindataNotMissing[nindex];
	           
		double[] gammaObservedSum_nindex = gammaObservedSum[nindex];

		for (int ns = 0; ns < numstates; ns++)
	        {
	           //going through each state
		   double[][] gammaksum_nseq_ns = gammaksum_nseq[ns];
	           double gammaObservedSum_nindex_ns = gammaObservedSum_nindex[ns];
		   for (int nmark = 0; nmark < numdatasets; nmark++)
	           {
		      //going through each mark
		      if (traindataNotMissing_nindex[nmark])
		      {
		         //only update non-missing
		         if (traindataObservedValues_nindex[nmark])
		         {
		            //updates the gamma sum for each mark when in state and observed 1
		            gammaksum_nseq_ns[nmark][1] += gammaObservedSum_nindex_ns;			    
			 }
		         else
		         {
		    	    //updates the gamma sum for each mark when in state and observed 0
		            gammaksum_nseq_ns[nmark][0] += gammaObservedSum_nindex_ns;			    
			 }
		      }
		   }
		}
	     }
	  }
	  //System.out.println("==>"+nseq+"\t"+dloglikeseq);
          dloglikeA[nseq] = dloglikeseq;	
	}
    }


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////





















    boolean[] threadslots;

    Object objlock = new Object();
    int nlaunched;


    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    /**
     * This is the core procedure for learning the parameters of the model loading from data
     */
    public void trainParametersParallelWithLoad() throws IOException
    {
        NumberFormat nf3 = NumberFormat.getInstance();
        nf3.setMaximumFractionDigits(3);
	nf3.setGroupingUsed(false);
	nf3.setMinimumFractionDigits(3);

        NumberFormat nf1 = NumberFormat.getInstance();
	nf1.setMaximumFractionDigits(1);
	nf1.setMinimumFractionDigits(1);
	nf1.setGroupingUsed(false);

       int niteration = 1;

       boolean bconverged = false;

       double dzerotransitioncutoff = Math.pow(10,-nzerotransitionpower);

       //number of non-zero transition for the 
       int nsparsecutoff = (int) (numstates * ChromHMM.SPARSECUTOFFRATIO);

       //number of non-zero transition that need to be less than this at the looser cut-off
       int nsparsecutofflooser = (int) (numstates * ChromHMM.SPARSECUTOFFLOOSERRATIO);

       double dprevloglike;

       //stores the maximum number of locations in any sequence and in each sequence
       //int[] numtime = new int[chromfiles.length];//traindataObservedIndex.length];
       int nmaxtime = 0;
       for (int nseq = 0; nseq < numtime.length; nseq++) //traindataObservedIndex.length;
       {
	   //numtime[nseq] = traindataObservedIndex[nseq].length;
          if (numtime[nseq] > nmaxtime)
	  {
      	     nmaxtime = numtime[nseq];
	  }
       }

       if (ChromHMM.BVERBOSE)
       {
          System.out.println("Maximum number of locations\t"+nmaxtime);
       }


       int numprocessors;
       //    for (int nseq = 0; nseq < traindataObservedIndex.length; nseq++)

       if (nmaxprocessors <= 0)
       {
	   numprocessors = Math.min(numtime.length,Runtime.getRuntime().availableProcessors()); //changed from  traindataObservedIndex.length;
       }
       else
       {
	   numprocessors = Math.min(numtime.length, //changed from  traindataObservedIndex.length;
				    Math.min(nmaxprocessors,Runtime.getRuntime().availableProcessors()));
       }

       //if (ChromHMM.BVERBOSE)
       {
	   if (numprocessors == 1)
              System.out.println("Using "+numprocessors+" thread for Baum-Welch training");
	   else
	       System.out.println("Using "+numprocessors+" threads for Baum-Welch training");
       }
       threadslots = new boolean[numprocessors];

       int ncurrnumincludeseq;
       int[] samples= null;
       if (numincludeseq >= 1)
       {
	   ncurrnumincludeseq = Math.min(numincludeseq, chromfiles.length);
	   samples = new int[ncurrnumincludeseq];
       }
       else
       {
	   ncurrnumincludeseq = chromfiles.length;
       }

       //stores the emission probability for the i^th combination of marks in the j^th state
       double[][][] emissionproducts_Pool = new double[numprocessors][nmaxtime][numstates]; //traindataObservedValues.length

       //stores the scale for the emission probability for the i^th combination of marks
       double[][] emissionproducts_scale_Pool;

       if (bscaleemissions)
       {
          emissionproducts_scale_Pool = new double[numprocessors][nmaxtime];
       }
       else
       {
          emissionproducts_scale_Pool = new double[numprocessors][];
       }


       //stores temporary product terms
       double[][] tempproductbetaemiss_Pool = new double[numprocessors][numstates];

       //This stores the alpha values at each time point and number of states
       double[][][] alpha_Pool = new double[numprocessors][nmaxtime][numstates];

       //Temporary storage of the gamma's for each state
       double[][] gamma_nt_Pool = new double[numprocessors][numstates];

       //Temporary storage of the beta values for each state
       double[][] beta_nt_Pool = new double[numprocessors][numstates];

       //Temporary storage of the beta values for each state at the next time point
       double[][] beta_ntp1_Pool = new double[numprocessors][numstates];

       //stores the scaling value for each time point
       double[][] scale_Pool = new double[numprocessors][nmaxtime];

       //stores the transition probabilities for each column
       double[][] coltransitionprobs = new double[numstates][numstates];

       //stores the sufficient statistic for the initital probability in each state for the last visit
       //double[][] gammainitstore = new double[chromfiles.length][numstates];
       double[][] gammainitstore = new double[ncurrnumincludeseq][numstates];

       //stores the sufficient statistics for computing the transition probabilities cumulated for each iteration
       double[][][] sxistore = new double[ncurrnumincludeseq][numstates][numstates];
       //double[][][] sxistore = new double[chromfiles.length][numstates][numstates];

       //stores the sufficient statistic for computing the emission probabilities
       double[][][][] gammaksumstore = 
	   //new double[chromfiles.length][numstates][numdatasets][numbuckets];
                 new double[ncurrnumincludeseq][numstates][numdatasets][numbuckets];

       //temporary storage in computation of sxi
       double[][][] sumforsxi_Pool = new double[numprocessors][numstates][numstates];

       //stores the sum of the gamma values associated with each combination in each state
       double[][][] gammaObservedSum_Pool = new double[numprocessors][nmaxtime][numstates]; //nmaxtime instead of traindataObservedValues.length

       //stores the indicies of the data
       int[][] traindataObservedIndex_Pool = new int[numprocessors][nmaxtime];

       //saving the mapping of signatures and chromsome observed on
       //stores whether there is a present call at each location
       boolean[][][] traindataObservedValues_Pool = new boolean[numprocessors][nmaxtime][numdatasets];

       //stores whether the mark is not considered missing
       boolean[][][] traindataNotMissing_Pool = new boolean[numprocessors][nmaxtime][numdatasets]; //usually nobserved


       double[] dloglikeA = new double[ncurrnumincludeseq];//chromfiles.length];

       int nelim = 0;
       long ltimeitr= System.currentTimeMillis();
       dprevloglike = Double.NEGATIVE_INFINITY;





       boolean[] bincludeseq = new boolean[chromfiles.length];
       for (int nk = 0; nk < bincludeseq.length; nk++)
       {
	   bincludeseq[nk] = true;
       }

       do
       {

	   if (numincludeseq >=1)
	   {
	       for (int ni = 0; ni < samples.length; ni++)
	       {
	          samples[ni] = ni;
	       }

	       //Random theRandom2 = new Random(412);//433

	       //sampling without replacement data locations wanted
	       for (int ni = ncurrnumincludeseq; ni < chromfiles.length; ni++)
	       {
	          if (theRandom.nextDouble() < ncurrnumincludeseq/((double) ni +1))
		  {
		     samples[theRandom.nextInt(ncurrnumincludeseq)] = ni;
		  }
	       }

               for (int nk = 0; nk < bincludeseq.length; nk++)
	       {
	           bincludeseq[nk] = false;
	       }

	       for (int nk = 0; nk < samples.length; nk++)
	       {
		   bincludeseq[samples[nk]] = true;
		   //System.out.println("including\t"+samples[nk]+"\t"+bincludeseq[samples[nk]]);
	       }
	   }

	   //dloglike= 0;	  

	   /*
	   //for (int ni = 0; ni < emissionproducts.length; ni++)
	  for (int ni = 0; ni < emissionproducts.length; ni++)
          {
             //going through each combination of marks
	     //this signature of marks is observed on the current chromosome so
	     //updating its emission probabilities
             double[] emissionproducts_ni = emissionproducts[ni];
	     boolean[] traindataObservedValues_ni = traindataObservedValues[ni];
	     boolean[] traindataNotMissing_ni = traindataNotMissing[ni];		  

	     boolean ballzero = true;

             for (int ns = 0; ns < numstates; ns++)
	     {
                double dproduct = 1;
	        double[][] emissionprobs_ns = emissionprobs[ns];

		for (int nmod = 0; nmod < numdatasets; nmod++)
	        {
	     	   if (traindataNotMissing_ni[nmod])
		   {
		      //we are include this marks emission probability
	              if (traindataObservedValues_ni[nmod])
		      {
		         //System.out.println("positive\t"+ns+"\t"+nmod+"\t1\t"+emissionprobs_ns[nmod][1]);
		         dproduct *= emissionprobs_ns[nmod][1];
		      }
		      else 
	              {
		         ///System.out.println("negative\t"+ns+"\t"+nmod+"\t0\t"+emissionprobs_ns[nmod][0]);
		         dproduct *= emissionprobs_ns[nmod][0];
		      }
		   }
		   // otherwise treated as missing omitting from product
		}
	        //System.out.println(ns+"\t"+dproduct);
	        emissionproducts_ni[ns] = dproduct;

		if (dproduct >= EPSILONEMISSIONS)
		{
	      	   ballzero = false;
	        }
	     }

	     if (ballzero)
	     {
	        for (int ns = 0; ns < numstates; ns++)
	        {
	       	   emissionproducts_ni[ns] = EPSILONEMISSIONS;
	        }
	     }
	  }
	   */

	  //stores in coltransitionprobs the transpose of transitionprobs
          for (int ni = 0; ni < numstates; ni++)
	  {
             double[] coltransitionprobs_ni = coltransitionprobs[ni];
	     for (int nj = 0; nj < numstates; nj++)
	     {
                coltransitionprobs_ni[nj] = transitionprobs[nj][ni];
             }
          }	     

          synchronized (objlock)
	  {
	      nlaunched = ncurrnumincludeseq;// chromfiles.length;// traindataObservedIndex.length;
	  }

          //for (int nseq = 0; nseq < traindataObservedIndex.length; nseq++)
	  int nincludeindex = 0;
          for (int nseq = 0; nseq < chromfiles.length; nseq++)
          {

	     if (bincludeseq[nseq])
	     {
		double[][] sxi_nseq = sxistore[nincludeindex];//nseq];
	        int numtime_nseq = numtime[nseq];
	        String chromfiles_nseq = chromfiles[nseq];
	        //int[] traindataObservedIndex_nseq = traindataObservedIndex[nseq];
	        //boolean[] traindataObservedSeqFlags_nseq = traindataObservedSeqFlags[nseq];
	        double[][][] gammaksum_nseq = gammaksumstore[nincludeindex];//nseq];
	        double[] gammainitstore_nseq = gammainitstore[nincludeindex];//nseq];
	       
	        NewThreadWithLoad myNewThreadWithLoad = new NewThreadWithLoad(
					 chromfiles_nseq,
					 traindataObservedIndex_Pool,
				      	 traindataObservedValues_Pool,
					 traindataNotMissing_Pool,
					 //traindataObservedIndex_nseq,
					 //traindataObservedSeqFlags_nseq,
					 gammaksum_nseq,
					 sxi_nseq,
					 numtime_nseq,
					 gammainitstore_nseq,
					 gammaObservedSum_Pool,
					 alpha_Pool,
					 emissionproducts_Pool,
					 gamma_nt_Pool,
					 coltransitionprobs,
					 scale_Pool,
					 beta_nt_Pool,
					 beta_ntp1_Pool,
					 tempproductbetaemiss_Pool,
					 sumforsxi_Pool,
					 nsparsecutoff,
					 nsparsecutofflooser,
					 dloglikeA,
					 nincludeindex,
                                         emissionproducts_scale_Pool);
		   nincludeindex++;

	           new Thread(myNewThreadWithLoad).start();
	     }	    
	  }	
	   

	  synchronized(objlock)
	  {
	     while (nlaunched > 0)
	     {
                try
                {
	           objlock.wait();
		}
	        catch (InterruptedException ex) {}
	     }
	  }


	  //normal EM
	  //executes the M-step after any pass through a sequence after one pass has been made through all sequences
          double dsum = 0;

	  //updating the initial probabilities
          for (int ni = 0; ni < numstates; ni++)
	  {
	     double dgammainitsum = 0;
             //for (int nitr = 0; nitr < traindataObservedIndex.length; nitr++)
             //for (int nitr = 0; nitr < chromfiles.length; nitr++)
             for (int nitr = 0; nitr < ncurrnumincludeseq; nitr++)
	     {
                dgammainitsum += gammainitstore[nitr][ni];
	     }

	     if ((bpseudo)&&(!bnopseudolast))
	     {
		dgammainitsum++;
	     }
	     probinit[ni] = dgammainitsum;
             dsum += dgammainitsum;
	  }
 
	  for (int ni = 0; ni < numstates; ni++)
          {
             probinit[ni] /= dsum;		
          }

	  //this indicates if there is a change on the set of 0 probability transitions
	  boolean bchange = false;
	  for (int ni = 0; ni < transitionprobs.length; ni++)
          {
	      dsum = 0;
	      //computes the denominator for the transition probabilities

	      int[] transitionprobsindex_ni = transitionprobsindex[ni];
	      double[] transitionprobs_ni =  transitionprobs[ni];
	      int transitionprobsnum_ni = transitionprobsnum[ni];
	      for (int nj = 0; nj < transitionprobsnum_ni; nj++)
              {
		  int ntransitionprobsindex_ni_nj = transitionprobsindex_ni[nj];
		  double dsxistoreitr = 0;
	          //for (int nitr = 0; nitr < chromfiles.length; nitr++) //traindataObservedIndex.length
	          for (int nitr = 0; nitr < ncurrnumincludeseq; nitr++) //traindataObservedIndex.length
	          {
		      //if (bincludeseq[nitr])
		     {
	                dsxistoreitr += sxistore[nitr][ni][ntransitionprobsindex_ni_nj];
		     }
	          }

		  if (bpseudo)
		  {
		      //pseudo-count to avoid divide by zero
		      dsxistoreitr++;
		  }
	          transitionprobs_ni[ntransitionprobsindex_ni_nj] = dsxistoreitr;

	          dsum += dsxistoreitr;
	      }
	       
              for (int nj = 0; nj < transitionprobsnum_ni; nj++)
	      {
		  int ntransitionprobsindex_ni_nj = transitionprobsindex_ni[nj];
		  //computes the updated transition probabilities

	          transitionprobs_ni[ntransitionprobsindex_ni_nj] /= dsum;
	    
		  if ((transitionprobs_ni[ntransitionprobsindex_ni_nj] < dzerotransitioncutoff) && (ni != ntransitionprobsindex_ni_nj))
	          {
		      //if falls below threshold eliminate the transition probabilities
		      elim[ni][ntransitionprobsindex_ni_nj] = true;
	              bchange = true;
	              nelim++;
	              transitionprobs_ni[ntransitionprobsindex_ni_nj] = 0;
	           }
	      }
	  }
	       
	      
          if (bchange)
	  {
	     //a transition was eliminated we need to update the probabilities
	     for (int ni = 0; ni < transitionprobs.length; ni++)
             {
		 int nindex = 0;		        
		 double ddenom = 0;
		 boolean[] elim_ni = elim[ni];
		 double[] transitionprobs_ni = transitionprobs[ni];
		 int[] transitionprobsindex_ni = transitionprobsindex[ni];
		 for (int nj = 0; nj < transitionprobs_ni.length; nj++)
	         {
	            if (!elim_ni[nj])
		    {
		       //we have not eliminated this transition
		       //stores its index in order and add sum to denominator
		       transitionprobsindex_ni[nindex] = nj;
	               ddenom += transitionprobs_ni[nj];
	               nindex++;
		    }
		 }

		 //renormalize the transition probabilities by the sum of the non-eliminated transitions
	         for (int nj = 0; nj < transitionprobs_ni.length; nj++)
	         {
	            transitionprobs_ni[nj] /= ddenom;
	         }
		 //update the number of valid transitions
	         transitionprobsnum[ni] = nindex; 
	     }

	     for (int ni = 0; ni < transitionprobs.length; ni++)
             {
		 int nindex =0;
		 int[] transitionprobsindexCol_ni = transitionprobsindexCol[ni];
		 for (int nj = 0; nj < transitionprobs[ni].length; nj++)
	         {
		    if (!elim[nj][ni])
	            {
		       //copy into the column of i the index of all non-eliminated transitions of i
	               transitionprobsindexCol_ni[nindex] = nj;
	               nindex++;
		    }
		 }
	         //updates the number of non-zero transitions from column i
	         transitionprobsnumCol[ni] = nindex; 
	     }
	  }
	    
	  //updating the emission parameters
	  for (int ns = 0; ns < numstates; ns++)
          {
	      double[][] emissionprobs_ns = emissionprobs[ns];


	      for (int nmark = 0; nmark < emissionprobs_ns.length; nmark++)
	      {
		  double[] emissionprobs_ns_nmark = emissionprobs_ns[nmark];
		  //can't used a general gamma sum because of missing emission vals
		  double dgammadenom = 0;

	          //updates gamma sum
                  for (int nbucket = 0; nbucket < numbuckets; nbucket++)
		  {
		      emissionprobs_ns_nmark[nbucket] = 0;
		      //ncurrnumincludeseq
	              //for (int nitr = 0; nitr < chromfiles.length; nitr++) //was traindataobserved length
		      //double dtempsum = 0;
	              for (int nitr = 0; nitr < ncurrnumincludeseq; nitr++) //was traindataobserved length
	              {
			  //if (bincludeseq[nitr])
			  //{
	                 emissionprobs_ns_nmark[nbucket] += gammaksumstore[nitr][ns][nmark][nbucket];
		          //}
	              }

		      if ((bpseudo) && ((!bnopseudolast)||(nmark<emissionprobs_ns.length-1)))
		      {
			  emissionprobs_ns_nmark[nbucket]++;
		      }
		      dgammadenom += emissionprobs_ns_nmark[nbucket];
		  }

		  //added to avoid NA
		  if (dgammadenom > 0)
	          {
	             for (int nbucket = 0; nbucket < numbuckets; nbucket++)
	             {
	                emissionprobs_ns_nmark[nbucket] /= dgammadenom;
		     }
		  }
	      }
	  }
       

          if (ChromHMM.BVERBOSE)
          {
             System.out.println("\t"+niteration+"\t"+dloglike);
	  }	     


	  dloglike = 0;
	  for (int nindex = 0; nindex < dloglikeA.length; nindex++)
	  {
	      //if (bincludeseq[nindex])
	      {
	         dloglike += dloglikeA[nindex];
	      }
	  }
    
	  double ddiff =(dloglike-dprevloglike);


	  //dconvergediff is only enforced if greater thanor equal to 0
          dprevloglike = dloglike;       

	  if (borderrows)
	  {
	     makeStateOrdering();
	  }

	  if (bordercols)
	  {
	     makeColOrdering();
	  }
	  //updates after each iteration the current status of the search
          printTransitionTable(niteration);
          printEmissionTable(niteration);

	  if (bprintimage)
	  {
             printEmissionImage(niteration);
             printTransitionImage(niteration);
	  }
	  printParametersToFile(niteration);

	  //we just completed a full iteration
          long ltimefinal =  System.currentTimeMillis();	  
	  double dtimechange = (ltimefinal-ltimeitr)/(double) 1000;
          bconverged = (((niteration >= nmaxiterations)||((ddiff< dconvergediff)&&(dconvergediff>=0)))||((dtimechange>nmaxseconds)&&(nmaxseconds>=0)));	  
          if (ChromHMM.BVERBOSE)
	  {
	     System.out.println(niteration+"\tTime Iteration\t"+dtimechange+"\t"+"\tElim\t"+nelim);
	     System.out.println("Full "+niteration+"\t"+dloglike+"\t"+dprevloglike+"\t"+ddiff);        
	  }


	  if (niteration == 1)
	  {
	      System.out.format("%10s %25s %10s %20s%n","Iteration","Estimated Log Likelihood", "Change","Total Time (secs)");
	      System.out.format("%10s %25s %10s %20s%n",""+niteration,""+nf3.format(dloglike),"-",""+nf1.format(dtimechange));
	  }
	  else
	  {
	      //System.out.format(niteration+"            "+nf3.format(dloglike)+"           "+nf3.format(ddiff)+"       "+nf1.format(dtimechange));
	      System.out.format("%10s %25s %10s %20s%n",""+niteration,""+nf3.format(dloglike),""+nf3.format(ddiff),""+nf1.format(dtimechange));
	  }
	  niteration++;
       }
       while (!bconverged);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    /**
     * This is the core procedure for learning the parameters of the model
     */
    public void trainParametersParallel() throws IOException
    {
        NumberFormat nf3 = NumberFormat.getInstance();
        nf3.setMaximumFractionDigits(3);
	nf3.setGroupingUsed(false);
	nf3.setMinimumFractionDigits(3);

        NumberFormat nf1 = NumberFormat.getInstance();
	nf1.setMaximumFractionDigits(1);
	nf1.setMinimumFractionDigits(1);
	nf1.setGroupingUsed(false);

       int niteration = 1;

       boolean bconverged = false;

       double dzerotransitioncutoff = Math.pow(10,-nzerotransitionpower);

       //number of non-zero transition for the 
       int nsparsecutoff = (int) (numstates * ChromHMM.SPARSECUTOFFRATIO);

       //number of non-zero transition that need to be less than this at the looser cut-off
       int nsparsecutofflooser = (int) (numstates * ChromHMM.SPARSECUTOFFLOOSERRATIO);

       double dprevloglike;

       //stores the maximum number of locations in any sequence and in each sequence
       int[] numtime = new int[traindataObservedIndex.length];
       int nmaxtime = 0;
       for (int nseq = 0; nseq < traindataObservedIndex.length; nseq++)
       {
          numtime[nseq] = traindataObservedIndex[nseq].length;
          if (numtime[nseq] > nmaxtime)
	  {
      	     nmaxtime = numtime[nseq];
	  }
       }

       if (ChromHMM.BVERBOSE)
       {
          System.out.println("Maximum number of locations\t"+nmaxtime);
       }


       int numprocessors;
       //    for (int nseq = 0; nseq < traindataObservedIndex.length; nseq++)

       if (nmaxprocessors <= 0)
       {
	   numprocessors = Math.min(traindataObservedIndex.length,Runtime.getRuntime().availableProcessors());
       }
       else
       {
	   numprocessors = Math.min(traindataObservedIndex.length,
				    Math.min(nmaxprocessors,Runtime.getRuntime().availableProcessors()));
       }

       //if (ChromHMM.BVERBOSE)
       {
	   if (numprocessors == 1)
              System.out.println("Using "+numprocessors+" thread for Baum-Welch training");
	   else
	       System.out.println("Using "+numprocessors+" threads for Baum-Welch training");
       }
       threadslots = new boolean[numprocessors];

       int ncurrnumincludeseq;
       int[] samples= null;
       if (numincludeseq >= 1)
       {
	   ncurrnumincludeseq = Math.min(numincludeseq, chromfiles.length);
	   samples = new int[ncurrnumincludeseq];
       }
       else
       {
	   ncurrnumincludeseq = chromfiles.length;
       }

       //stores the emission probability for the i^th combination of marks in the j^th state
       double[][] emissionproducts = new double[traindataObservedValues.length][numstates];


       //stores the scale for the emission probability for the i^th combination of marks
       double[] emissionproducts_scale = null;

       if (bscaleemissions)
       {
          emissionproducts_scale = new double[traindataObservedValues.length];
       }

       //stores temporary product terms
       double[][] tempproductbetaemiss_Pool = new double[numprocessors][numstates];

       //This stores the alpha values at each time point and number of states
       double[][][] alpha_Pool = new double[numprocessors][nmaxtime][numstates];

       //Temporary storage of the gamma's for each state
       double[][] gamma_nt_Pool = new double[numprocessors][numstates];

       //Temporary storage of the beta values for each state
       double[][] beta_nt_Pool = new double[numprocessors][numstates];

       //Temporary storage of the beta values for each state at the next time point
       double[][] beta_ntp1_Pool = new double[numprocessors][numstates];

       //stores the scaling value for each time point
       double[][] scale_Pool = new double[numprocessors][nmaxtime];

       //stores the transition probabilities for each column
       double[][] coltransitionprobs = new double[numstates][numstates];

       //stores the sufficient statistic for the initital probability in each state for the last visit
       double[][] gammainitstore = //new double[traindataObservedIndex.length][numstates];
 	              new double[ncurrnumincludeseq][numstates];

       //stores the sufficient statistics for computing the transition probabilities cumulated for each iteration
       double[][][] sxistore = //new double[traindataObservedIndex.length][numstates][numstates];
 	             new double[ncurrnumincludeseq][numstates][numstates];

       //stores the sufficient statistic for computing the emission probabilities
       double[][][][] gammaksumstore = 
	         new double[ncurrnumincludeseq][numstates][numdatasets][numbuckets];
                  //new double[traindataObservedIndex.length][numstates][numdatasets][numbuckets];

       //temporary storage in computation of sxi
       double[][][] sumforsxi_Pool = new double[numprocessors][numstates][numstates];

       //stores the sum of the gamma values associated with each combination in each state
       double[][][] gammaObservedSum_Pool = new double[numprocessors][traindataObservedValues.length][numstates];


       double[] dloglikeA = new double[ncurrnumincludeseq]; 
	   //new double[traindataObservedIndex.length];

       int nelim = 0;
       long ltimeitr= System.currentTimeMillis();
       dprevloglike = Double.NEGATIVE_INFINITY;


       boolean[] bincludeseq = new boolean[chromfiles.length];
       for (int nk = 0; nk < bincludeseq.length; nk++)
       {
	   bincludeseq[nk] = true;
       }


       do
       {
	   //dloglike= 0;	  

	   if (numincludeseq >=1)
	   {
	       for (int ni = 0; ni < samples.length; ni++)
	       {
	          samples[ni] = ni;
	       }

	       //Random theRandom2 = new Random(412);//433

	       //sampling without replacement data locations wanted
	       for (int ni = ncurrnumincludeseq; ni < chromfiles.length; ni++)
	       {
	          if (theRandom.nextDouble() < ncurrnumincludeseq/((double) ni +1))
		  {
		     samples[theRandom.nextInt(ncurrnumincludeseq)] = ni;
		  }
	       }

               for (int nk = 0; nk < bincludeseq.length; nk++)
	       {
	           bincludeseq[nk] = false;
	       }

	       for (int nk = 0; nk < samples.length; nk++)
	       {
		   bincludeseq[samples[nk]] = true;
		   //System.out.println("including\t"+samples[nk]+"\t"+bincludeseq[samples[nk]]);
	       }
	   }


           if (bscaleemissions)
	   {
	      for (int ni = 0; ni < emissionproducts.length; ni++)
	      {
	         //going through each combination of marks
		 //this signature of marks is observed on the current chromosome so
	         //updating its emission probabilities
		 double[] emissionproducts_ni = emissionproducts[ni];
	         boolean[] traindataObservedValues_ni = traindataObservedValues[ni];
		 boolean[] traindataNotMissing_ni = traindataNotMissing[ni];

		 for (int ns = 0; ns < numstates; ns++)
		 {
	            emissionproducts_ni[ns] = 1;
		 }
		 emissionproducts_scale[ni] = 0;

		 for (int nmod = 0; nmod < numdatasets; nmod++)
		 {
		    for (int ns = 0; ns < numstates; ns++)
		    {
		       if (traindataNotMissing_ni[nmod])
		       {
		          //we are include this marks emission probability
			  if (traindataObservedValues_ni[nmod])
			  {
			     emissionproducts_ni[ns] *= emissionprobs[ns][nmod][1];
			  }
			  else
			  {
			     emissionproducts_ni[ns] *= emissionprobs[ns][nmod][0];
			  }
		       }
		       // otherwise treated as missing omitting from product
		    }

		    double dmaxval = 0;
		    for (int ns = 0; ns < numstates; ns++)
		    {
		       if (emissionproducts_ni[ns] > dmaxval)
		       {
			   dmaxval = emissionproducts_ni[ns];
		       }
		    }

		    if (dmaxval <= 0)//EPSILONEMISSIONS)
		    {
		       for (int ns = 0; ns < numstates; ns++)
		       {
		          emissionproducts_ni[ns] = 1;
		       }
		       //approximating log-likelihood
		       emissionproducts_scale[ni] += Math.log(EPSILONEMISSIONS);
		    }
		    else
		    {
		       for (int ns = 0; ns < numstates; ns++)
		       {
			  emissionproducts_ni[ns]/= dmaxval;
		       }
		       emissionproducts_scale[ni] += Math.log(dmaxval);
		    }
		 }
	      }
	   }
	   else
	   {
	      for (int ni = 0; ni < emissionproducts.length; ni++)
              {
                 //going through each combination of marks
	         //this signature of marks is observed on the current chromosome so
	         //updating its emission probabilities
                 double[] emissionproducts_ni = emissionproducts[ni];
	         boolean[] traindataObservedValues_ni = traindataObservedValues[ni];
	         boolean[] traindataNotMissing_ni = traindataNotMissing[ni];		  

	         boolean ballzero = true;

                 for (int ns = 0; ns < numstates; ns++)
	         {
                    double dproduct = 1;
	            double[][] emissionprobs_ns = emissionprobs[ns];

		    for (int nmod = 0; nmod < numdatasets; nmod++)
	            {
	     	       if (traindataNotMissing_ni[nmod])
		       {
		          //we are include this marks emission probability
	                  if (traindataObservedValues_ni[nmod])
		          {
		             //System.out.println("positive\t"+ns+"\t"+nmod+"\t1\t"+emissionprobs_ns[nmod][1]);
		             dproduct *= emissionprobs_ns[nmod][1];
		          }
		          else 
	                  {
		             ///System.out.println("negative\t"+ns+"\t"+nmod+"\t0\t"+emissionprobs_ns[nmod][0]);
		             dproduct *= emissionprobs_ns[nmod][0];
			  }
		       }
		       // otherwise treated as missing omitting from product
		    }
	            //System.out.println(ns+"\t"+dproduct);
	            emissionproducts_ni[ns] = dproduct;

		    if (dproduct >= EPSILONEMISSIONS)
		    {
	      	       ballzero = false;
		    }
		 }

	         if (ballzero)
	         {
	            for (int ns = 0; ns < numstates; ns++)
	            {
	       	       emissionproducts_ni[ns] = EPSILONEMISSIONS;
		    }
		 }
	      }
	   }


	  //stores in coltransitionprobs the transpose of transitionprobs
          for (int ni = 0; ni < numstates; ni++)
	  {
             double[] coltransitionprobs_ni = coltransitionprobs[ni];
	     for (int nj = 0; nj < numstates; nj++)
	     {
                coltransitionprobs_ni[nj] = transitionprobs[nj][ni];
             }
          }	     

          synchronized (objlock)
	  {
	      nlaunched = ncurrnumincludeseq;//traindataObservedIndex.length;
	  }


	  int nincludeindex = 0;
          for (int nseq = 0; nseq < traindataObservedIndex.length; nseq++)
          {
	      //System.out.println(nseq+"\t"+bincludeseq[nseq]);
	     if (bincludeseq[nseq])
	     {
		double[][] sxi_nseq = sxistore[nincludeindex];//[nseq]; //bug fixed here in v1.15 when using -n option with lowmem
	        int numtime_nseq = numtime[nseq];
	        int[] traindataObservedIndex_nseq = traindataObservedIndex[nseq];
	        boolean[] traindataObservedSeqFlags_nseq = traindataObservedSeqFlags[nseq];
	        double[][][] gammaksum_nseq = gammaksumstore[nincludeindex];// gammaksumstore[nseq];
	        double[] gammainitstore_nseq = gammainitstore[nincludeindex];// gammainitstore[nseq];

	        NewThread myNewThread = new NewThread(traindataObservedIndex_nseq,
					 traindataObservedSeqFlags_nseq,
					 gammaksum_nseq,
					 sxi_nseq,
					 numtime_nseq,
					 gammainitstore_nseq,
					 gammaObservedSum_Pool,
					 alpha_Pool,
					 emissionproducts,
					 gamma_nt_Pool,
					 coltransitionprobs,
					 scale_Pool,
					 beta_nt_Pool,
					 beta_ntp1_Pool,
					 tempproductbetaemiss_Pool,
					 sumforsxi_Pool,
					 nsparsecutoff,
					 nsparsecutofflooser,
					 dloglikeA,
				   	 nincludeindex,
                                         emissionproducts_scale);
						      //nseq);
		nincludeindex++;
	        new Thread(myNewThread).start();
	     } 
	  }	
	   

	  synchronized(objlock)
	  {
	     while (nlaunched > 0)
	     {
                try
                {
	           objlock.wait();
		}
	        catch (InterruptedException ex) {}
	     }
	  }

	  //normal EM
	  //executes the M-step after any pass through a sequence after one pass has been made through all sequences
          double dsum = 0;

	  //updating the inital probabilities
          for (int ni = 0; ni < numstates; ni++)
	  {
	     double dgammainitsum = 0;
             //for (int nitr = 0; nitr < traindataObservedIndex.length; nitr++)
             for (int nitr = 0; nitr < ncurrnumincludeseq; nitr++)
	     {
		 //if (bincludeseq[nitr])
		{
                   dgammainitsum += gammainitstore[nitr][ni];
		}
	     }

	     if ((bpseudo)&&(!bnopseudolast))
	     {
		 dgammainitsum++;
	     }
	     probinit[ni] = dgammainitsum;
             dsum += dgammainitsum;
	  }
 
	  for (int ni = 0; ni < numstates; ni++)
          {
             probinit[ni] /= dsum;		
          }

	  //this indicates if there is a change on the set of 0 probability transitions
	  boolean bchange = false;
	  for (int ni = 0; ni < transitionprobs.length; ni++)
          {
	      dsum = 0;
	      //computes the denominator for the transition probabilities

	      int[] transitionprobsindex_ni = transitionprobsindex[ni];
	      double[] transitionprobs_ni =  transitionprobs[ni];
	      int transitionprobsnum_ni = transitionprobsnum[ni];
	      for (int nj = 0; nj < transitionprobsnum_ni; nj++)
              {
		  int ntransitionprobsindex_ni_nj = transitionprobsindex_ni[nj];
		  double dsxistoreitr = 0;
	          //for (int nitr = 0; nitr < traindataObservedIndex.length; nitr++)
	          for (int nitr = 0; nitr < ncurrnumincludeseq; nitr++)
	          {
		      //if (bincludeseq[nitr])
		     {
	                dsxistoreitr += sxistore[nitr][ni][ntransitionprobsindex_ni_nj];
		     }
	          }

		  if (bpseudo)
		  {
		      //pseudo-count to avoid divide by zero
		      dsxistoreitr++;
		  }
	          transitionprobs_ni[ntransitionprobsindex_ni_nj] = dsxistoreitr;

	          dsum += dsxistoreitr;
	      }
	       
              for (int nj = 0; nj < transitionprobsnum_ni; nj++)
	      {
		  int ntransitionprobsindex_ni_nj = transitionprobsindex_ni[nj];
		  //computes the updated transition probabilities
	          transitionprobs_ni[ntransitionprobsindex_ni_nj] /= dsum;
	    
		  if ((transitionprobs_ni[ntransitionprobsindex_ni_nj] < dzerotransitioncutoff) && (ni != ntransitionprobsindex_ni_nj))
	          {
		      //if falls below threshold eliminate the transition probabilities
		      elim[ni][ntransitionprobsindex_ni_nj] = true;
	              bchange = true;
	              nelim++;
	              transitionprobs_ni[ntransitionprobsindex_ni_nj] = 0;
	           }
	      }
	  }
	       
	      
          if (bchange)
	  {
	     //a transition was eliminated we need to update the probabilities
	     for (int ni = 0; ni < transitionprobs.length; ni++)
             {
		 int nindex = 0;		        
		 double ddenom = 0;
		 boolean[] elim_ni = elim[ni];
		 double[] transitionprobs_ni = transitionprobs[ni];
		 int[] transitionprobsindex_ni = transitionprobsindex[ni];
		 for (int nj = 0; nj < transitionprobs_ni.length; nj++)
	         {
	            if (!elim_ni[nj])
		    {
		       //we have not eliminated this transition
		       //stores its index in order and add sum to denominator
		       transitionprobsindex_ni[nindex] = nj;
	               ddenom += transitionprobs_ni[nj];
	               nindex++;
		    }
		 }

		 //renormalize the transition probabilities by the sum of the non-eliminated transitions
	         for (int nj = 0; nj < transitionprobs_ni.length; nj++)
	         {
	            transitionprobs_ni[nj] /= ddenom;
	         }
		 //update the number of valid transitions
	         transitionprobsnum[ni] = nindex; 
	     }

	     for (int ni = 0; ni < transitionprobs.length; ni++)
             {
		 int nindex =0;
		 int[] transitionprobsindexCol_ni = transitionprobsindexCol[ni];
		 for (int nj = 0; nj < transitionprobs[ni].length; nj++)
	         {
		    if (!elim[nj][ni])
	            {
		       //copy into the column of i the index of all non-eliminated transitions of i
	               transitionprobsindexCol_ni[nindex] = nj;
	               nindex++;
		    }
		 }
	         //updates the number of non-zero transitions from column i
	         transitionprobsnumCol[ni] = nindex; 
	     }
	  }
	    
	  //updating the emission parameters
	  for (int ns = 0; ns < numstates; ns++)
          {
	      double[][] emissionprobs_ns = emissionprobs[ns];

	      for (int nmark = 0; nmark < emissionprobs_ns.length; nmark++)
	      {
		  double[] emissionprobs_ns_nmark = emissionprobs_ns[nmark];
		  //can't used a general gamma sum because of missing emission vals
		  double dgammadenom = 0;

	          //updates gamma sum
                  for (int nbucket = 0; nbucket < numbuckets; nbucket++)
		  {
		      emissionprobs_ns_nmark[nbucket] = 0;
	                
	              //for (int nitr = 0; nitr < traindataObservedIndex.length; nitr++)
	              for (int nitr = 0; nitr < ncurrnumincludeseq; nitr++)
	              {
			  //if (bincludeseq[nitr])
			 {
	                    emissionprobs_ns_nmark[nbucket] += gammaksumstore[nitr][ns][nmark][nbucket];
			 }
	              }

		      if ((bpseudo) && ((!bnopseudolast)||(nmark<emissionprobs_ns.length-1)))
		      {
			  emissionprobs_ns_nmark[nbucket]++;
		      }
		      dgammadenom += emissionprobs_ns_nmark[nbucket];
		  }

		  //added to avoid NA
		  if (dgammadenom > 0)
		  {
	             for (int nbucket = 0; nbucket < numbuckets; nbucket++)
	             {
	                emissionprobs_ns_nmark[nbucket] /= dgammadenom;
	             }
		  }
	      }
	  }
       

          if (ChromHMM.BVERBOSE)
          {
             System.out.println("\t"+niteration+"\t"+dloglike);
	  }	     


	  dloglike = 0;
	  for (int nindex = 0; nindex < dloglikeA.length; nindex++)
	  {
	      //if (bincludeseq[nindex])
	      {
	         dloglike += dloglikeA[nindex];
	      }
	  }
    
	  double ddiff =(dloglike-dprevloglike);


	  //dconvergediff is only enforced if greater thanor equal to 0
          dprevloglike = dloglike;       

	  if (borderrows)
	  {
	     makeStateOrdering();
	  }

	  if (bordercols)
	  {
	     makeColOrdering();
	  }
	  //updates after each iteration the current status of the search
          printTransitionTable(niteration);
          printEmissionTable(niteration);
	  if (bprintimage)
	  {
             printEmissionImage(niteration);
             printTransitionImage(niteration);
	  }
	  printParametersToFile(niteration);

	  //we just completed a full iteration
          long ltimefinal =  System.currentTimeMillis();	  
	  double dtimechange = (ltimefinal-ltimeitr)/(double) 1000;
          bconverged = (((niteration >= nmaxiterations)||((ddiff< dconvergediff)&&(dconvergediff>=0)))||((dtimechange>nmaxseconds)&&(nmaxseconds>=0)));	  
          if (ChromHMM.BVERBOSE)
	  {
	     System.out.println(niteration+"\tTime Iteration\t"+dtimechange+"\t"+"\tElim\t"+nelim);
	     System.out.println("Full "+niteration+"\t"+dloglike+"\t"+dprevloglike+"\t"+ddiff);        
	  }


	  if (niteration == 1)
	  {
	      System.out.format("%10s %25s %10s %20s%n","Iteration","Estimated Log Likelihood", "Change","Total Time (secs)");
	      System.out.format("%10s %25s %10s %20s%n",""+niteration,""+nf3.format(dloglike),"-",""+nf1.format(dtimechange));
	  }
	  else
	  {
	      //System.out.format(niteration+"            "+nf3.format(dloglike)+"           "+nf3.format(ddiff)+"       "+nf1.format(dtimechange));
	      System.out.format("%10s %25s %10s %20s%n",""+niteration,""+nf3.format(dloglike),""+nf3.format(ddiff),""+nf1.format(dtimechange));
	  }
	  niteration++;
       }
       while (!bconverged);
    }


    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    /**
     * Loads in the input data
     * If there are multiple cell type associated with the files that should be
     * indicated by a prefix before an '_'
     */
    public void loadDataFileStubs() throws IOException
    {

	if (szinputfilelist == null)
        {
	    //takes all files in the directory with a _binary
	   File dir = new File(szinputdir);
           String[] chromfilesall = dir.list();
	   if (chromfilesall == null)
	   {
	       throw new IllegalArgumentException(szinputdir+" is not a valid directory!");
	   }
	   ArrayList alfiles = new ArrayList();
	   for (int nfile = 0; nfile < chromfilesall.length; nfile++)
	   {
	       if ((chromfilesall[nfile].contains("_binary"))&&(!(new File(chromfilesall[nfile])).isHidden()))
	       {
	           //added hidden check in v.1.11
	          alfiles.add(chromfilesall[nfile]);
	       }
	   }

	   if (alfiles.size() == 0)
	   {
	       throw new IllegalArgumentException("No files found in "+szinputdir+" containing '_binary' that are not Hidden");
	   }

 	   //stores them in chromfiles
           chromfiles = new String[alfiles.size()];
           for (int nfile = 0; nfile < chromfiles.length; nfile++)
           {
              chromfiles[nfile] = (String) alfiles.get(nfile);
           }	   
	}
	else
	{
	    //loads in the input coords list
	    BufferedReader brfiles = Util.getBufferedReader(szinputfilelist);

	    ArrayList alfiles = new ArrayList();
	    String szLine;
	    while ((szLine = brfiles.readLine())!=null)
	    {
		alfiles.add(szLine);
	    }
	    brfiles.close(); 

	    //stores them in chromfiles
	    chromfiles = new String[alfiles.size()];
	    for (int nfile = 0; nfile < chromfiles.length; nfile++)
	    {
		chromfiles[nfile] = (String) alfiles.get(nfile);
	    }
	}

	Arrays.sort(chromfiles);//gives a deterministic reproducible starting order to the chromfiles
	
	//randomly orders the chromosome files to visit
	RecIntDouble[] recA = new RecIntDouble[chromfiles.length];
	String[] tempchromfiles = new String[chromfiles.length];

	//now going to randomize the order of chromosomes if theRandom is not null
	//the order of chromosome matters when doing an incremental expectation maximization
	for (int ni = 0; ni < chromfiles.length; ni++)
	{
	   tempchromfiles[ni] = chromfiles[ni];
	   if (theRandom == null)
	   {
	       recA[ni] = new RecIntDouble(ni, ni);
	   }
	   else
	   {
	      recA[ni] = new RecIntDouble(ni,theRandom.nextDouble());
	   }
	}

	if (theRandom != null)
	{
	    //already in order if random is null
	    Arrays.sort(recA,new RecIntDoubleCompare());
	}

	cellSeq = new String[chromfiles.length];
	chromSeq = new String[chromfiles.length];

	for (int ni = 0; ni < chromfiles.length; ni++)
	{
	    //swapping into chromfiles the sorted chromosome ordering
	    chromfiles[ni] = tempchromfiles[recA[ni].nindex];
	}
	
	traindataObservedIndex = new int[chromfiles.length][]; //number of columns depends on number of lines in file
        numtime = new int[chromfiles.length];

	//HashMap hmObserved = new HashMap(); //maps an observation string to an index and set of flags
 
	int nobserved = 0;
        //PrintWriter pw = null;

	for (int nfile = 0; nfile < chromfiles.length; nfile++)
        {
	    if (ChromHMM.BVERBOSE)
	    {
	       System.out.println("reading\t"+szinputdir+" "+chromfiles[nfile]);
	    }
	    BufferedReader br = Util.getBufferedReader(szinputdir+"/"+chromfiles[nfile]);
	    String szLine = br.readLine(); //first line tells cell type and chromosome
	    if (szLine == null)
	    {
		throw new IllegalArgumentException(szinputdir+"/"+chromfiles[nfile]+" is empty!");
	    }
	    StringTokenizer st = new StringTokenizer(szLine,"\t");
            if (!st.hasMoreTokens())
	    {
	       throw new IllegalArgumentException("First line must contain cell type and chromosome. No entries found.");
	    }

	    cellSeq[nfile] = st.nextToken().trim();//added trim in version 1.23
	    if (!st.hasMoreTokens())
	    {
		throw new IllegalArgumentException("First line must contain cell type and chromosome. Only one entry found.");
	    }
	    chromSeq[nfile] = st.nextToken().trim();//added trim in version 1.23

	    if (st.hasMoreTokens())
	    {
		throw new IllegalArgumentException("First line should only contain cell type and chromosome");
	    }
	    szLine = br.readLine(); //reading header
            //to output binary
	    if (szLine == null)
	    {
		throw new IllegalArgumentException(szinputdir+"/"+chromfiles[nfile]+" only has one line!");
	    }
	    st = new StringTokenizer(szLine,"\t");
	    int numtokens = st.countTokens();
	    if (nfile == 0)
	    {
		//first time reading header taking tokens
	       datasets = new String[numtokens];
	       int ntoken = 0;
	       while (st.hasMoreTokens())
	       {
		  datasets[ntoken] = st.nextToken().trim();
	          ntoken++;
	       }
 	       //numdatasets is the number of marks we are integrating
      	       numdatasets = datasets.length; 
	    }
	    else
	    {
		//Requires number of tokens to match
		if (numtokens != datasets.length)
		{
		    throw new IllegalArgumentException(" found a file with header with "+numtokens+" entries, which does not match another with "+datasets.length); //updated in v1.20 to be more informative
		}

		int ntoken = 0;
		//Gives warning if a header column does not match
	        while (st.hasMoreTokens())
	        {
		   String sztoken = st.nextToken().trim();
	           if (!datasets[ntoken].equals(sztoken))
		   {
		       System.out.println("WARNING headers do not match between "+chromfiles[nfile]+" and "+chromfiles[0]);
		   }
		   ntoken++;
	        }
	    }



	    // LineNumberReader  lnr = new LineNumberReader(br);//new FileReader(new File("File1")));
	    //while (lnr.skip(Long.MAX_VALUE) > 0);
	    //int nlinecount = lnr.getLineNumber();


	    
	    int nlinecount  = 0;
	    while ((szLine = br.readLine())!=null)
            {
	    	nlinecount++;
	    }
	    
	    //System.out.println("num lines is\t"+nlinecount+"\t"+nfile);
	    br.close();
	    //lnr.close();
	    numtime[nfile] = nlinecount;


	}
    }


    /*
	    ArrayList aldata = new ArrayList();
	    while ((szLine = br.readLine())!=null)
	    {
		st = new StringTokenizer(szLine,"\t");
		StringBuffer sb = new StringBuffer();
		
		for (int ncol = 0; ncol < numdatasets; ncol++)
		{

		    if (!st.hasMoreTokens())
		    {
			throw new IllegalArgumentException("Found line without "+numdatasets+" values in file "+chromfiles[nfile]);
		    }

		    String sztoken = st.nextToken();
		    
		    if (sztoken.equals("0"))
		    {
			sb.append("0");
		    }
		    else if (sztoken.equals("1"))
		    {
			sb.append("1");
		    }
		    else if (sztoken.equals("2"))
		    {
			//this means missing
			sb.append("2");
		    }
		    else
		    {
			throw new IllegalArgumentException("Unrecognized value "+sztoken+" found in "+szinputdir+"/"+chromfiles[nfile]);
		    }
		}

		aldata.add(sb.toString());
	    }
	    br.close();

	    int nsize = aldata.size();
	    traindataObservedIndex[nfile] = new int[nsize];
	    int[] traindataObservedIndex_nfile = traindataObservedIndex[nfile];

	    for (int nrow = 0; nrow < nsize; nrow++)
	    {
		BigInteger theBigInteger = new BigInteger((String) aldata.get(nrow),3);
		ObservedRec theObservedRec  = (ObservedRec) hmObserved.get(theBigInteger);
		boolean[] flagA;

		if (theObservedRec == null)
		{
		    //this is the first time we encountered this combination of marks
		    flagA = new boolean[chromfiles.length];
		    //recording which chromsomes this mark combination was observed
		    flagA[nfile] =true;

		    //System.out.println(szmappingbyte.length());
		    //storing a mapping from observed byte string to an integer index in alFlags and alObserved
		    hmObserved.put(theBigInteger, new ObservedRec(nobserved,flagA));

		    //saving this observed index
		    traindataObservedIndex[nfile][nrow] = nobserved;

		    //increments the number of observed combinations of marks
		    nobserved++;
		}
		else
		{
		    //updating that this signature was observed on this chromosome
		    theObservedRec.flagA[nfile] = true;
		    //storing the index of the flags associated with this row 
		    traindataObservedIndex_nfile[nrow] = theObservedRec.nobserved;
		}
	    }
	}
	    
	//saving the mapping of signatures and chromsome observed on

	//stores whether there is a present call at each location
	traindataObservedValues = new boolean[nobserved][numdatasets];

	//stores whether the mark is not considered missing
	traindataNotMissing = new boolean[nobserved][numdatasets];

	//stores whether this sequence combination appears on the chromosome
	traindataObservedSeqFlags = new boolean[chromfiles.length][traindataObservedValues.length];

	Iterator hmObservedIterator = hmObserved.entrySet().iterator();
	while (hmObservedIterator.hasNext())
	{
	    Map.Entry pairs = (Map.Entry) hmObservedIterator.next();
	    BigInteger theBigInteger = (BigInteger) pairs.getKey();
	    String szmapping = theBigInteger.toString(3);  //getting back the mapping string

	    ObservedRec theObservedRec = (ObservedRec) pairs.getValue();
	    int ncurrindex = theObservedRec.nobserved;//this is an index on which obervation combination it is

	    boolean[] traindataObservedValues_ncurrindex = traindataObservedValues[ncurrindex];
	    boolean[] traindataNotMissing_ncurrindex = traindataNotMissing[ncurrindex]; 
	    
	    //if the mapping string is less than the number of data sets then 
	    //there are leading 0's will set for leading 0's not missing and absent
	    int numch = szmapping.length();
	    int numleading0 = numdatasets - numch;
	    for (int nj = 0; nj < numleading0; nj++)
	    {
		traindataObservedValues_ncurrindex[nj] = false;
		traindataNotMissing_ncurrindex[nj] = true;
	    }

	    int nmappedindex = numleading0; //starting from the leading 0 position
	    for (int nj = 0; nj < numch; nj++)
	    {
	       char ch = szmapping.charAt(nj);

	       if (ch == '0')
	       {
		   traindataObservedValues_ncurrindex[nmappedindex] = false;
		   traindataNotMissing_ncurrindex[nmappedindex] = true;
	       }
	       else if (ch=='1')
	       {
		   traindataObservedValues_ncurrindex[nmappedindex] = true;
		   traindataNotMissing_ncurrindex[nmappedindex] = true;
	       }
	       else
	       {
		   //missing data
		   traindataObservedValues_ncurrindex[nmappedindex] = false;
		   traindataNotMissing_ncurrindex[nmappedindex] = false;
	       }
	       nmappedindex++;
	    }

	    boolean[] currFlags = theObservedRec.flagA;
	    for (int nj = 0; nj < chromfiles.length; nj++)
            {
		//storing at this observation whether it is found for each chromosome
		traindataObservedSeqFlags[nj][ncurrindex] = currFlags[nj];
	    }
	}	
    */




    ///////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * Loads in the input data
     * If there are multiple cell type associated with the files that should be 
     * indicated by a prefix before an '_'
     */
    public void loadData() throws IOException
    {	
	if (szinputfilelist == null)
        {
	    //takes all files in the directory with a _binary
	   File dir = new File(szinputdir);
           String[] chromfilesall = dir.list();
	   if (chromfilesall == null)
	   {
	       throw new IllegalArgumentException(szinputdir+" is not a valid directory!");
	   }
	   ArrayList alfiles = new ArrayList();
	   for (int nfile = 0; nfile < chromfilesall.length; nfile++)
	   {
	       if ((chromfilesall[nfile].contains("_binary"))&&(!(new File(chromfilesall[nfile])).isHidden()))
	       {
	           //added hidden check in v.1.11
	          alfiles.add(chromfilesall[nfile]);
	       }
	   }

	   if (alfiles.size() == 0)
	   {
	       throw new IllegalArgumentException("No files found in "+szinputdir+" containing '_binary' that are not Hidden");
	   }

 	   //stores them in chromfiles
           chromfiles = new String[alfiles.size()];
           for (int nfile = 0; nfile < chromfiles.length; nfile++)
           {
              chromfiles[nfile] = (String) alfiles.get(nfile);
           }	   
	}
	else
	{
	    //loads in the input coords list
	    BufferedReader brfiles = Util.getBufferedReader(szinputfilelist);

	    ArrayList alfiles = new ArrayList();
	    String szLine;
	    while ((szLine = brfiles.readLine())!=null)
	    {
		alfiles.add(szLine);
	    }
	    brfiles.close(); 

	    //stores them in chromfiles
	    chromfiles = new String[alfiles.size()];
	    for (int nfile = 0; nfile < chromfiles.length; nfile++)
	    {
		chromfiles[nfile] = (String) alfiles.get(nfile);
	    }
	}

	Arrays.sort(chromfiles);//gives a deterministic reproducible starting order to the chromfiles       
	
	//randomly orders the chromosome files to visit
	RecIntDouble[] recA = new RecIntDouble[chromfiles.length];
	String[] tempchromfiles = new String[chromfiles.length];

	//now going to randomize the order of chromosomes if theRandom is not null
	//the order of chromosome matters when doing an incremental expectation maximization
	for (int ni = 0; ni < chromfiles.length; ni++)
	{
	   tempchromfiles[ni] = chromfiles[ni];
	   if (theRandom == null)
	   {
	       recA[ni] = new RecIntDouble(ni, ni);
	   }
	   else
	   {
	      recA[ni] = new RecIntDouble(ni,theRandom.nextDouble());
	   }
	}

	if (theRandom != null)
	{
	    //already in order if random is null
	    Arrays.sort(recA,new RecIntDoubleCompare());
	}

	cellSeq = new String[chromfiles.length];
	chromSeq = new String[chromfiles.length];

	for (int ni = 0; ni < chromfiles.length; ni++)
	{
	    //swapping into chromfiles the sorted chromosome ordering
	    chromfiles[ni] = tempchromfiles[recA[ni].nindex];
	}
	
	traindataObservedIndex = new int[chromfiles.length][]; //number of columns depends on number of lines in file

	HashMap hmObserved = new HashMap(); //maps an observation string to an index and set of flags
 
	int nobserved = 0;
        //PrintWriter pw = null;

	for (int nfile = 0; nfile < chromfiles.length; nfile++)
        {
	    if (ChromHMM.BVERBOSE)
	    {
	       System.out.println("reading\t"+szinputdir+" "+chromfiles[nfile]);
	    }
	    BufferedReader br = Util.getBufferedReader(szinputdir+"/"+chromfiles[nfile]);
	    String szLine = br.readLine(); //first line tells cell type and chromosome
	    if (szLine == null)
	    {
		throw new IllegalArgumentException(szinputdir+"/"+chromfiles[nfile]+" is empty!");
	    }
	    StringTokenizer st = new StringTokenizer(szLine,"\t");
            if (!st.hasMoreTokens())
	    {
	       throw new IllegalArgumentException("First line must contain cell type and chromosome. No entries found.");
	    }

	    cellSeq[nfile] = st.nextToken().trim(); //added trim in version 1.23
	    if (!st.hasMoreTokens())
	    {
		throw new IllegalArgumentException("First line must contain cell type and chromosome. Only one entry found.");
	    }
	    chromSeq[nfile] = st.nextToken().trim(); //added trim in version 1.23

	    if (st.hasMoreTokens())
	    {
		throw new IllegalArgumentException("First line should only contain cell type and chromosome");
	    }
	    szLine = br.readLine(); //reading header
            //to output binary
	    if (szLine == null)
	    {
		throw new IllegalArgumentException(szinputdir+"/"+chromfiles[nfile]+" only has one line!");
	    }
	    st = new StringTokenizer(szLine,"\t");
	    int numtokens = st.countTokens();
	    if (nfile == 0)
	    {
		//first time reading header taking tokens
	       datasets = new String[numtokens];
	       int ntoken = 0;
	       while (st.hasMoreTokens())
	       {
		  datasets[ntoken] = st.nextToken().trim();
		  ntoken++;
	       }
	    }
	    else
	    {
		//Requires number of tokens to match
		if (numtokens != datasets.length)
		{
		    throw new IllegalArgumentException(" found a file with header with "+numtokens+" entries, which does not match another with "+datasets.length); //updated in v1.20 to be more informative
		}

		int ntoken = 0;
		//Gives warning if a header column does not match
	        while (st.hasMoreTokens())
	        {
		   String sztoken = st.nextToken().trim();
	           if (!datasets[ntoken].equals(sztoken))
		   {
		       System.out.println("WARNING headers do not match between "+chromfiles[nfile]+" and "+chromfiles[0]);
		   }
		   ntoken++;
	        }
	    }

	    //numdatasets is the number of marks we are integrating
      	    numdatasets = datasets.length; 

	    ArrayList aldata = new ArrayList();
	    while ((szLine = br.readLine())!=null)
	    {
		st = new StringTokenizer(szLine,"\t ");
		StringBuffer sb = new StringBuffer();
		
		for (int ncol = 0; ncol < numdatasets; ncol++)
		{

		    if (!st.hasMoreTokens())
		    {
			throw new IllegalArgumentException("Found line without "+numdatasets+" values in file "+chromfiles[nfile]);
		    }

		    String sztoken = st.nextToken();//.trim();
		    
		    if (sztoken.equals("0"))
		    {
			sb.append("0");
		    }
		    else if (sztoken.equals("1"))
		    {
			sb.append("1");
		    }
		    else if (sztoken.equals("2"))
		    {
			//this means missing
			sb.append("2");
		    }
		    else
		    {
			throw new IllegalArgumentException("Unrecognized value "+sztoken+" found in "+szinputdir+"/"+chromfiles[nfile]);
		    }
		}

		aldata.add(sb.toString());
	    }
	    br.close();

	    int nsize = aldata.size();
	    traindataObservedIndex[nfile] = new int[nsize];
	    int[] traindataObservedIndex_nfile = traindataObservedIndex[nfile];

	    for (int nrow = 0; nrow < nsize; nrow++)
	    {
		BigInteger theBigInteger = new BigInteger((String) aldata.get(nrow),3);
		ObservedRec theObservedRec  = (ObservedRec) hmObserved.get(theBigInteger);
		boolean[] flagA;

		if (theObservedRec == null)
		{
		    //this is the first time we encountered this combination of marks
		    flagA = new boolean[chromfiles.length];
		    //recording which chromsomes this mark combination was observed
		    flagA[nfile] =true;

		    //System.out.println(szmappingbyte.length());
		    //storing a mapping from observed byte string to an integer index in alFlags and alObserved
		    hmObserved.put(theBigInteger, new ObservedRec(nobserved,flagA));

		    //saving this observed index
		    traindataObservedIndex_nfile[nrow] = nobserved; //change in 1.20 for efficiency to be_nfile

		    //increments the number of observed combinations of marks
		    nobserved++;
		}
		else
		{
		    //updating that this signature was observed on this chromosome
		    theObservedRec.flagA[nfile] = true;
		    //storing the index of the flags associated with this row 
		    traindataObservedIndex_nfile[nrow] = theObservedRec.nobserved;
		}
	    }
	}
	    
	//saving the mapping of signatures and chromsome observed on

	//stores whether there is a present call at each location
	traindataObservedValues = new boolean[nobserved][numdatasets];

	//stores whether the mark is not considered missing
	traindataNotMissing = new boolean[nobserved][numdatasets];

	//stores whether this sequence combination appears on the chromosome
	traindataObservedSeqFlags = new boolean[chromfiles.length][traindataObservedValues.length];

	Iterator hmObservedIterator = hmObserved.entrySet().iterator();
	while (hmObservedIterator.hasNext())
	{
	    Map.Entry pairs = (Map.Entry) hmObservedIterator.next();
	    BigInteger theBigInteger = (BigInteger) pairs.getKey();
	    String szmapping = theBigInteger.toString(3);  //getting back the mapping string

	    ObservedRec theObservedRec = (ObservedRec) pairs.getValue();
	    int ncurrindex = theObservedRec.nobserved;//this is an index on which obervation combination it is

	    boolean[] traindataObservedValues_ncurrindex = traindataObservedValues[ncurrindex];
	    boolean[] traindataNotMissing_ncurrindex = traindataNotMissing[ncurrindex]; 
	    
	    //if the mapping string is less than the number of data sets then 
	    //there are leading 0's will set for leading 0's not missing and absent
	    int numch = szmapping.length();
	    int numleading0 = numdatasets - numch;
	    for (int nj = 0; nj < numleading0; nj++)
	    {
		traindataObservedValues_ncurrindex[nj] = false;
		traindataNotMissing_ncurrindex[nj] = true;
	    }

	    int nmappedindex = numleading0; //starting from the leading 0 position
	    for (int nj = 0; nj < numch; nj++)
	    {
	       char ch = szmapping.charAt(nj);

	       if (ch == '0')
	       {
		   traindataObservedValues_ncurrindex[nmappedindex] = false;
		   traindataNotMissing_ncurrindex[nmappedindex] = true;
	       }
	       else if (ch=='1')
	       {
		   traindataObservedValues_ncurrindex[nmappedindex] = true;
		   traindataNotMissing_ncurrindex[nmappedindex] = true;
	       }
	       else
	       {
		   //missing data
		   traindataObservedValues_ncurrindex[nmappedindex] = false;
		   traindataNotMissing_ncurrindex[nmappedindex] = false;
	       }
	       nmappedindex++;
	    }

	    boolean[] currFlags = theObservedRec.flagA;
	    for (int nj = 0; nj < chromfiles.length; nj++)
            {
		//storing at this observation whether it is found for each chromosome
		traindataObservedSeqFlags[nj][ncurrindex] = currFlags[nj];
	    }
	}	
    }


    //////////////////////////////////////////////////////////////////////////////////////////////
  
    public static void main(String[] args) throws IOException
    {

	boolean bok = true;
	String szcommand = "";
	if (args.length >= 1)
	{
	   szcommand = args[0];
	}
	else
	{
	    bok = false;
	}

	if (szcommand.equalsIgnoreCase("Version"))
	{
	    System.out.println("This is Version 1.26 of ChromHMM (c) Copyright 2008-2012 Massachusetts Institute of Technology");
	}
        else if ((szcommand.equals("BinarizeBam"))||(szcommand.equalsIgnoreCase("BinarizeBed")))
	{
	    boolean bstacked = false;
	    boolean bgzip = false;
	    boolean bpairend = false;
	    boolean bmixed = false;
	    boolean bsplit = false;
	    String szcontroldir=null;
	    int nflankwidthcontrol = 5;
	    int nshift = 100;
	    boolean bcenterinterval = false;
	    int noffsetleft = 0;
	    int noffsetright = 1;
	    double dpoissonthresh = 0.0001;
	    double dfoldthresh = 0;
	    double dcountthresh = 0;
	    boolean bcontainsthresh = true;
	    int npseudocount = 1;
	    int nbinsize = ChromHMM.DEFAULT_BINSIZEBASEPAIRS;
	    String szcolfields=null;
	    String szoutputcontroldir=null;
	    String szoutputsignaldir = null;
	    int npseudocountcontrol = 1;
	    boolean bpeaks = false;
	    int numsplitbins = ChromHMM.DEFAULT_NUMSPLITBINS;
	    int nsplitindex = -1;
	    boolean bshift = false;
	    boolean bflagsplitbins = false;

	    boolean bflagsplitcols = false;

	    boolean bsplitcols = false;
	    int nnummarksplit = 10; 
	    int nmarksplitindex = -1; 

	    int nargindex = 1;
	    if (args.length <= 4)
	    {
		bok = false;
	    }
	    else
	    {
               try
	       {
	          while (nargindex < args.length-4)
	          {
                     if (args[nargindex].equals("-b"))
		     {
		        nbinsize = Integer.parseInt(args[++nargindex]);
		     }
		     else if (args[nargindex].equals("-c"))
		     {
		        szcontroldir = args[++nargindex];
		     }
		     else if ((args[nargindex].equals("-colfields"))&&(szcommand.equalsIgnoreCase("BinarizeBed")))
		     {
		        szcolfields = args[++nargindex];
		     }
		     else if (args[nargindex].equals("-e"))
		     {
		        noffsetright= Integer.parseInt(args[++nargindex]);
		     }
		     else if (args[nargindex].equals("-f"))
		     {
		        dfoldthresh = Double.parseDouble(args[++nargindex]);
		     }
                     else if (args[nargindex].equals("-g"))
		     {
		        dcountthresh = Double.parseDouble(args[++nargindex]);
		     }
                     else if (args[nargindex].equals("-gzip"))
		     {
		         bgzip = true;
		     }
                     else if (args[nargindex].equals("-i"))
		     {
			 nsplitindex = Integer.parseInt(args[++nargindex]);
		     }
		     else if (args[nargindex].equals("-j"))
		     {
                         numsplitbins = Integer.parseInt(args[++nargindex]);
			 bflagsplitbins = true;
			 //bsplit = true;
		     }
		     else if (args[nargindex].equals("-k"))
		     {
                         nmarksplitindex = Integer.parseInt(args[++nargindex]);
			 bflagsplitcols = true;
		     }
		     else if (args[nargindex].equals("-m"))
		     {
			 nnummarksplit = Integer.parseInt(args[++nargindex]);
			 bflagsplitcols = true;
		     }
		     else if (args[nargindex].equals("-n"))
		     {
		         nshift = Integer.parseInt(args[++nargindex]);
			 bshift = true;
		     }
  	  	     else if (args[nargindex].equals("-o"))
		     {
		        szoutputcontroldir = args[++nargindex];
			File f = new File(szoutputcontroldir);
			if (!f.exists())
			{
			    if (!f.mkdirs())
			    {
				throw new IllegalArgumentException(szoutputcontroldir+" does not exist and could not be created!");
			    }
			}
		     }
		     else if (args[nargindex].equals("-p"))
		     {
		        dpoissonthresh = Double.parseDouble(args[++nargindex]);
		     }
		     else if (args[nargindex].equals("-s"))
		     {
		        noffsetleft = Integer.parseInt(args[++nargindex]);
		     }
		     else if (args[nargindex].equals("-splitrows"))
		     {
			 bsplit = true;
		     }
                     else if (args[nargindex].equals("-splitcols"))
		     {
		         bsplitcols = true;
		     }
                     else if (args[nargindex].equals("-stacked"))
		     {
		        bstacked = true;
		     }
		     else if (args[nargindex].equals("-strictthresh"))
		     {
		        bcontainsthresh = false;
		     }
		     else if (args[nargindex].equals("-t"))
		     {
		        szoutputsignaldir = args[++nargindex];
			File f = new File(szoutputsignaldir);
			if (!f.exists())
			{
			    if (!f.mkdirs())
			    {
				throw new IllegalArgumentException(szoutputsignaldir+" does not exist and could not be created!");
			    }
			}
		     }
		     else if (args[nargindex].equals("-u"))
		     {
		        npseudocountcontrol = Integer.parseInt(args[++nargindex]);
		     }
		     else if (args[nargindex].equals("-center"))
		     {
		        bcenterinterval = true;
		     }
		     else if (args[nargindex].equals("-peaks"))
		     {
			 bpeaks = true;
		     }
		     else if ((args[nargindex].equals("-paired"))&&(szcommand.equalsIgnoreCase("BinarizeBam")))
		     {
			 bpairend = true;
		     }
		     else if ((args[nargindex].equals("-mixed"))&&(szcommand.equalsIgnoreCase("BinarizeBam")))
		     {
			 bmixed = true;
		     }
		     else if (args[nargindex].equals("-w"))
		     {
		        nflankwidthcontrol = Integer.parseInt(args[++nargindex]);
		     }
		     else
		     {
		        bok = false;
		        break;
		     }
		     nargindex++;
		  }
	       }
	       catch (NumberFormatException ex)
	       {
	          bok = false;
	       } 	       
	    }

	    if ((bpairend) && (bmixed||bcenterinterval||bshift||bpeaks))
	    {
		bok = false;
	    }

	    if ((!bsplit)&& (bflagsplitbins))
	    {
	       bok = false;
	    }

	    if ((!bsplitcols)&&(bflagsplitcols))
	    {
		bok = false;
	    }



	    if ((bok)&&(nargindex == args.length-4))
	    {
	       String szchromlengthfile = args[nargindex++];
	       String szmarkdir = args[nargindex++];
	       String szcellmarkfiletable = args[nargindex++];
	       String szoutputbinarydir= args[nargindex];

	       File f = new File(szoutputbinarydir);
	       if (!f.exists())
	       {
	          if (!f.mkdirs())
	          {
		     throw new IllegalArgumentException(szoutputbinarydir+" does not exist and could not be created!");
		  }
	       }

	       if (szcontroldir == null)
	       {
		   szcontroldir = szmarkdir;
	       }



	       int nmaxsplit = 0;
	       if (bsplitcols)
	       {
	          String szLineCellMark;
		  int ncellmarkline = 0;

		  BufferedReader brcellmark = Util.getBufferedReader(szcellmarkfiletable);
		  HashSet hscombo = new HashSet();

		  while ((szLineCellMark=brcellmark.readLine()) != null)
	          {
		     if (szLineCellMark.trim().equals("")) continue;

		     StringTokenizer stcellmark = new StringTokenizer(szLineCellMark,"\t");

		     if (stcellmark.countTokens() < 3)
		     {
		        throw new IllegalArgumentException("In "+szcellmarkfiletable+" "+szLineCellMark+" does not have at least three columns");
		     }
		     String szcell = stcellmark.nextToken().trim(); //added trim in v1.20 to remove leading and trailing white space
		     String szmark = stcellmark.nextToken().trim();

		     if (hscombo.contains(szcell+"\t"+szmark))
		     {
	                throw new IllegalArgumentException("In "+szcellmarkfiletable+" "+szcell+"\t"+szmark+" found twice, but -splitcols specified");
		     }
		     hscombo.add(szcell+"\t"+szmark);

		     if (ncellmarkline % nnummarksplit == 0)
		     {
	                nmaxsplit++;
		     }
		     ncellmarkline++;
		  }
		  brcellmark.close();
	       }

	       if (nsplitindex >= 0)
	       {
	          if (bpeaks)
		  {
	             if ((bsplitcols)&&(nmarksplitindex==-1))
		     {
	                //need to figure out maxsplit

			for (nmarksplitindex = 0; nmarksplitindex < nmaxsplit; nmarksplitindex++)
		        {
			    String szoutputbinarydirsub = szoutputbinarydir+"/SET"+(nmarksplitindex+1);
			    f = new File(szoutputbinarydirsub);
			    if (!f.exists())
		       	    {
		      	       if (!f.mkdirs())
			       {
			          throw new IllegalArgumentException(szoutputbinarydirsub+" does not exist and could not be created!");
			       }
			    }

			    Preprocessing.makeBinaryDataFromPeaksSplit(szchromlengthfile, szmarkdir, szoutputbinarydir, szcellmarkfiletable,
					                         nbinsize, bgzip, numsplitbins, nsplitindex, noffsetleft,noffsetright,
								       bsplitcols,nnummarksplit,nmarksplitindex,bstacked);
			}
		     }
		     else
		     {
		        if (bsplitcols)
		        {
			   if (nmarksplitindex >= nmaxsplit)
			   {
			       throw new IllegalArgumentException("nmarksplitindex"+" value of "+nmarksplitindex+" is greater than the maximum of "+(nmaxsplit-1));
			   }
		           szoutputbinarydir = szoutputbinarydir+"/SET"+(nmarksplitindex+1);
		           f = new File(szoutputbinarydir);
		           if (!f.exists())
		           {
		              if (!f.mkdirs())
		              {
		                 throw new IllegalArgumentException(szoutputbinarydir+" does not exist and could not be created!");
			      }
			   }
			}
		     
		        Preprocessing.makeBinaryDataFromPeaksSplit(szchromlengthfile, szmarkdir, szoutputbinarydir, szcellmarkfiletable,
					                         nbinsize, bgzip, numsplitbins, nsplitindex, noffsetleft,noffsetright,
								   bsplitcols,nnummarksplit,nmarksplitindex, bstacked);
		     }
		  }
		  else
		  {
		     bok = false;
		  }
	       }
	       else
	       {

		  if ((bsplitcols)&&(nmarksplitindex==-1))
		  {
		      //going through all
		      //need to figure out maxsplit

		      for (nmarksplitindex = 0; nmarksplitindex < nmaxsplit; nmarksplitindex++)
		      {
			 String szoutputbinarydirsub = szoutputbinarydir+"/SET"+(nmarksplitindex+1);
		         f = new File(szoutputbinarydirsub);
		         if (!f.exists())
		         {
		            if (!f.mkdirs())
		            {
		               throw new IllegalArgumentException(szoutputbinarydirsub+" does not exist and could not be created!");
			    }
			 }

			 String szoutputsignaldirsub = null;
			 if (szoutputsignaldir != null)
			 {
			    szoutputsignaldirsub = szoutputsignaldir+"/SET"+(nmarksplitindex+1);
		            f = new File(szoutputsignaldirsub);
		            if (!f.exists())
		            {
		               if (!f.mkdirs())
		               {
		                  throw new IllegalArgumentException(szoutputsignaldirsub+" does not exist and could not be created!");
			       }
			    }
			 }

                         String szoutputcontroldirsub = null;
			 if (szoutputcontroldirsub != null)
			 {
			    szoutputcontroldirsub = szoutputcontroldir+"/SET"+(nmarksplitindex+1);
		            f = new File(szoutputcontroldirsub);
		            if (!f.exists())
		            {
		               if (!f.mkdirs())
		               {
		                  throw new IllegalArgumentException(szoutputcontroldirsub+" does not exist and could not be created!");
			       }
			    }
			 }
		  
	                 Preprocessing.makeBinaryDataFromBed(szchromlengthfile,szmarkdir,szcontroldir,nflankwidthcontrol,szcellmarkfiletable,
						nshift,bcenterinterval, noffsetleft,noffsetright,szoutputsignaldir,
						szoutputbinarydirsub,szoutputcontroldir,
					        dpoissonthresh,dfoldthresh,bcontainsthresh,
					        npseudocountcontrol,nbinsize,szcolfields,bpeaks, dcountthresh,szcommand.equalsIgnoreCase("BinarizeBam"),
							     bpairend, bgzip, bsplit, numsplitbins, bsplitcols,nnummarksplit,nmarksplitindex, bstacked, bmixed);	   	        
		      }
		  }
		  else
		  {
		     if (bsplitcols)
		     {
			if (nmarksplitindex >= nmaxsplit)
			{
			    throw new IllegalArgumentException("nmarksplitindex"+" value of "+nmarksplitindex+" is greater than the maximum of "+(nmaxsplit-1));
			}
		        szoutputbinarydir = szoutputbinarydir+"/SET"+(nmarksplitindex+1);
		        f = new File(szoutputbinarydir);
		        if (!f.exists())
		        {
		           if (!f.mkdirs())
		           {
		              throw new IllegalArgumentException(szoutputbinarydir+" does not exist and could not be created!");
			   }
			}

			if (szoutputsignaldir != null)
			{
			   szoutputsignaldir = szoutputsignaldir+"/SET"+(nmarksplitindex+1);
		           f = new File(szoutputsignaldir);
		           if (!f.exists())
		           {
		              if (!f.mkdirs())
		              { 
		                 throw new IllegalArgumentException(szoutputsignaldir+" does not exist and could not be created!");
			      }
			   }
			}

			if (szoutputcontroldir != null)
			{
			   szoutputcontroldir = szoutputcontroldir+"/SET"+(nmarksplitindex+1);
		           f = new File(szoutputcontroldir);
		           if (!f.exists())
		           {
		              if (!f.mkdirs())
		              {
		                 throw new IllegalArgumentException(szoutputcontroldir+" does not exist and could not be created!");
			      }
			   }
			}
		     }
		  
	             Preprocessing.makeBinaryDataFromBed(szchromlengthfile,szmarkdir,szcontroldir,nflankwidthcontrol,szcellmarkfiletable,
						nshift,bcenterinterval, noffsetleft,noffsetright,szoutputsignaldir,
						szoutputbinarydir,szoutputcontroldir,
					        dpoissonthresh,dfoldthresh,bcontainsthresh,
					        npseudocountcontrol,nbinsize,szcolfields,bpeaks, dcountthresh,szcommand.equalsIgnoreCase("BinarizeBam"),
							 bpairend, bgzip, bsplit, numsplitbins, bsplitcols,nnummarksplit,nmarksplitindex,bstacked,bmixed);	   	          
		   
		  }
	       }
	    }
	    else
            {
       	       bok = false;
	    }

	
	    if (!bok)
	    {
	       if (szcommand.equalsIgnoreCase("BinarizeBed"))
	       {
		   //v1.18 update
                  System.out.println("usage BinarizeBed [-b binsize][-c controldir][-center][-colfields chromosome,start,end[,strand]][-e offsetend][-f foldthresh]"+
                                  "[-g signalthresh][-gzip][-n shift][-o outputcontroldir][-p poissonthresh][-peaks [-i splitrowindex]][-s offsetstart][-splitcols [-k splitcolindex][-m numsplitcols]][-splitrows [-j numsplitbins]][-stacked][-strictthresh][-t outputsignaldir]"+
                                  "[-u pseudocountcontrol][-w flankwidthcontrol] "+
                                  "chromosomelengthfile inputbeddir cellmarkfiletable outputbinarydir");
	       }
	       else
	       {
		   System.out.println("usage BinarizeBam [-b binsize][-c controldir][-e offsetend][-f foldthresh]"+
                                  "[-g signalthresh][-gzip][[-o outputcontroldir][-p poissonthresh][-paired|[-mixed][-center][-n shift][-peaks [-i splitindex]]"+
                                  "[-s offsetstart][-splitcols [-k splitcolindex][-m numsplitcols]][-splitrows [-j numsplitbins]][-stacked][-strictthresh][-t outputsignaldir]"+
                                  "[-u pseudocountcontrol][-w flankwidthcontrol] "+
				      "chromosomelengthfile inputbamdir cellmarkfiletable outputbinarydir");
	       }
	    }
	}
	else if (szcommand.equalsIgnoreCase("BinarizeSignal"))
	{
	    boolean bcontainsthresh = true;
	    boolean bgzip = false;
	    boolean bsplit = false;
	    double dfoldthresh = 0;
	    double dcountthresh = 0;
	    double dpoissonthresh = 0.0001;
	    int nflankwidthcontrol = 5;
	    int npseudocountcontrol = 1;
	    int numsplitbins = ChromHMM.DEFAULT_NUMSPLITBINS;
	    boolean bflagsplitbins = false;

	    String szsignaldir = null;
	    String szoutputdir = null;
	    String szcontroldir = null;

	    if (args.length <= 2)
	    {
		bok = false;
	    }
	    else
	    {
	       int nargindex = 1;
               try
	       {
	          while (nargindex < args.length-2)
	          {
		     if (args[nargindex].equals("-c"))
		     {
		        szcontroldir = args[++nargindex];
		     }
		     else if (args[nargindex].equals("-f"))
		     {
		        dfoldthresh = Double.parseDouble(args[++nargindex]);
		     }
                     else if (args[nargindex].equals("-g"))
		     {
		        dcountthresh = Double.parseDouble(args[++nargindex]);
		     }
                     else if (args[nargindex].equals("-gzip"))
		     {
		        bgzip = true;
		     }
		     else if (args[nargindex].equals("-j"))
		     {
			 numsplitbins = Integer.parseInt(args[++nargindex]);
			 bflagsplitbins = true;
		     }
		     else if (args[nargindex].equals("-p"))
		     {
		        dpoissonthresh = Double.parseDouble(args[++nargindex]);
		     }
		     else if (args[nargindex].equals("-splitrows"))
		     {
			 bsplit = true;
		     }
		     else if (args[nargindex].equals("-strictthresh"))
		     {
		        bcontainsthresh = false;
		     }
		     else if (args[nargindex].equals("-u"))
		     {
		        npseudocountcontrol = Integer.parseInt(args[++nargindex]);
		     }
		     else if (args[nargindex].equals("-w"))
		     {
		        nflankwidthcontrol = Integer.parseInt(args[++nargindex]);
		     } 
		     else
		     {
		        bok = false;
			break;
		     }
		     nargindex++;
		  }
	       }
	       catch (NumberFormatException ex)
	       {
		   bok = false;
	       }

	       if ((!bsplit)&& (bflagsplitbins))
	       {
	          bok = false;
	       }
	    
	       if ((bok)&&(nargindex == args.length-2))
	       {
	          szsignaldir = args[nargindex++];
	          szoutputdir = args[nargindex];

		  File f = new File(szoutputdir);
	          if (!f.exists())
	          {
	             if (!f.mkdirs())
	             {
		        throw new IllegalArgumentException(szoutputdir+" does not exist and could not be created!");
		     }
		  }		  
	      
	          if (szcontroldir!=null)
	          {
                     Preprocessing.makeBinaryDataFromSignalAgainstControl(szsignaldir,szcontroldir, szoutputdir,
								   dpoissonthresh, dfoldthresh,bcontainsthresh, 
									  nflankwidthcontrol,npseudocountcontrol, dcountthresh,bgzip, bsplit, numsplitbins);
	          }
	          else
	          {
		     Preprocessing.makeBinaryDataFromSignalUniform(szsignaldir, szoutputdir, dpoissonthresh, 
								   dfoldthresh, bcontainsthresh,dcountthresh,bgzip, bsplit, numsplitbins);
	          }
	       }
	       else
	       {
		   bok = false;
	       }
	    }

	    if (!bok)
            {
               System.out.println("usage BinarizeSignal [-c controldir][-f foldthresh][-g signalthresh][-gzip][-p poissonthresh][-splitrows [-j numsplitbins]][-strictthresh][-u pseudocountcontrol][-w flankwidth] signaldir outputdir");
            }	       
	    
	}
	else if (szcommand.equalsIgnoreCase("MergeBinary"))
	{

	    boolean bgzip = false;
	    int numsplitbins = ChromHMM.DEFAULT_NUMSPLITBINS;
	    boolean bsplit = false;
	    String szdirlistfile = null;
	    String sztype = "binary";
	    String szinputdir = null;
	    String szoutputdir = null;
	    boolean bflagsplitbins = false;

	    if (args.length <= 2)
	    {
		bok = false;
	    }
	    else
	    {
	       int nargindex = 1;
               try
	       {
	          while (nargindex < args.length-2)
	          {
		     if (args[nargindex].equals("-f"))
		     {
			 szdirlistfile = args[++nargindex];
		     }
		     else if (args[nargindex].equals("-t"))
		     {
			 sztype = args[++nargindex];
		     }
		     else if (args[nargindex].equals("-j"))
		     {
			numsplitbins = Integer.parseInt(args[++nargindex]);
			bflagsplitbins = true;
		        //bsplit = true;
		     }
		     else if (args[nargindex].equals("-splitrows"))
		     {
			bsplit = true;
		     }
                     else if (args[nargindex].equals("-gzip"))
		     {
		        bgzip = true;
		     }
		     else
		     {
		        bok = false;
			break;
		     }
		     nargindex++;
		  }
	       }
	       catch (NumberFormatException ex)
	       {
		   bok = false;
	       }

	       if ((!bsplit)&& (bflagsplitbins))
	       {
	          bok = false;
	       }
	    
	       if ((bok)&&(nargindex == args.length-2))
	       {
	          szinputdir = args[nargindex++];
	          szoutputdir = args[nargindex];

		  File f = new File(szoutputdir);
	          if (!f.exists())
	          {
	             if (!f.mkdirs())
	             {
		        throw new IllegalArgumentException(szoutputdir+" does not exist and could not be created!");
		     }
		  }		  
		  Preprocessing.mergeBinarizedFiles(szinputdir, szoutputdir, szdirlistfile, bsplit, numsplitbins, bgzip, sztype);
	       }
	       else
	       {
		   bok = false;
	       }
	    }

	    if (!bok)
            {
               System.out.println("usage MergeBinary [-f dirlistfile][-gzip][-splitrows [-j numsplitbins]][-t type] inputdir outputdir");
            }	       
	}
	else if (szcommand.equalsIgnoreCase("CompareModels"))
	{
	    String szmainmodel;
	    String szinputdir;
	    String szoutputprefix;

	    int nr=ChromHMM.DEFAULTCOLOR_R;
	    int ng=ChromHMM.DEFAULTCOLOR_G;
	    int nb=ChromHMM.DEFAULTCOLOR_B;
	    boolean bprintimage = true;

	    int nargindex = 1;

	    if (args.length <= 2)
	    {
		bok = false;
	    }
	    else
	    {
               try
	       {
	          while (nargindex < args.length-3)
		  {
		     if (args[nargindex].equals("-color"))
		     {
		        String szcolor = args[++nargindex];
		        StringTokenizer stcolor = new StringTokenizer(szcolor,",");
		        if (stcolor.countTokens()==3)
		        {
		           nr = Integer.parseInt(stcolor.nextToken().trim());
			   ng = Integer.parseInt(stcolor.nextToken().trim());
		           nb = Integer.parseInt(stcolor.nextToken().trim());
			}
		        else
		        {
		           bok = false;
		        }
		     }
		     else if (args[nargindex].equals("-noimage"))
		     {
	                bprintimage = false;
		     }	
		     else
		     {
	                bok = false;
		     }
		     nargindex++;
		  }
	       }
	       catch (NumberFormatException ex)
	       {
	          bok = false;
	       }
	    }
	    
	    if (nargindex != args.length-3)
	    {
		bok = false;
	    }

	    if (bok)
	    {
	       szmainmodel = args[nargindex++];
	       szinputdir = args[nargindex++];
	       szoutputprefix = args[nargindex];
	       StateAnalysis.makeModelEmissionCompare(szmainmodel,szinputdir,szoutputprefix,new Color(nr,ng,nb),bprintimage);
	    }
	    
	    if (!bok)
	    {
	       System.out.println("usage CompareModels [-color r,g,b][-noimage] referencemodelemissions comparedir outputprefix");
	    }	    
	}
	else if (szcommand.equalsIgnoreCase("StatePruning"))
	{
	    boolean beuclidean = true;
	    int nindex = 1;
	    String szinputdir;
	    String szoutputdir;
	    if (nindex >= args.length)
	    {
		bok = false;
	    }
	    else 
	    {
	       if (args[nindex].equals("-correlation"))
	       {
		  beuclidean = false;
		  nindex++;
	       }
	       if (nindex+1 == (args.length-1))
	       {
	          szinputdir = args[nindex++];
	          szoutputdir = args[nindex];	
		  File f = new File(szoutputdir);
	          if (!f.exists())
	          {
	             if (!f.mkdirs())
	             {
		        throw new IllegalArgumentException(szoutputdir+" does not exist and could not be created!");
		     }
		  }

		  NestedEliminateInitialize.nestedEliminateInitialize(szinputdir, szoutputdir,beuclidean);
	       }
	       else
	       {
		   bok = false;
	       }
	    }
	
	    if (!bok)
            {
               System.out.println("usage: StatePruning [-correlation] inputdir outputdir"); 
            }	       	    
	}
	else if (szcommand.equalsIgnoreCase("EvalSubset"))
	{
	    int nr=ChromHMM.DEFAULTCOLOR_R;
	    int ng=ChromHMM.DEFAULTCOLOR_G;
	    int nb=ChromHMM.DEFAULTCOLOR_B;
	    Color theColor = new Color(nr, ng, nb);

	    String szinputfilelist = null;
	    boolean breadposterior = false;
	    boolean breadstatebyline = false;
	    boolean breadsegment = false;
	    boolean bprintimage = true;
            boolean bscaleemissions = false;
	    boolean bscalebeta = false;
	    String szchromlengthfile = null;
	    int nbinsize = ChromHMM.DEFAULT_BINSIZEBASEPAIRS;
	    String szoutfileID = "";
	    boolean bappend = false;
	    int nargindex = 1;
	    boolean blowmem = false;

            try
	    {
	       while (nargindex < args.length-5)
	       {
                  if (args[nargindex].equals("-color"))
		  {
		     String szcolor = args[++nargindex];
		     StringTokenizer stcolor = new StringTokenizer(szcolor,",");
		     if (stcolor.countTokens()==3)
		     {
		        nr = Integer.parseInt(stcolor.nextToken().trim());
		        ng = Integer.parseInt(stcolor.nextToken().trim());
		        nb = Integer.parseInt(stcolor.nextToken().trim());
		     }
		     else
		     {
		        bok = false;
		     }
		  }
	          else if (args[nargindex].equals("-b"))
		  {
		     nbinsize = Integer.parseInt(args[++nargindex]);
		  }
		  else if (args[nargindex].equals("-f"))
		  {
		     szinputfilelist = args[++nargindex];
		  }
		  else if (args[nargindex].equals("-i"))
		  {
		     szoutfileID = args[++nargindex];
		  }
		  else if (args[nargindex].equals("-lowmem"))
		  {
		      blowmem = true;
		  }
                  else if (args[nargindex].equals("-append"))
		  {
		      bappend = true;
		  }
                  else if (args[nargindex].equals("-many"))
	          {
		      bscaleemissions = true;
		  }
                  else if (args[nargindex].equals("-noimage"))
		  {
		      bprintimage = false;
		  }
		  else if (args[nargindex].equals("-readposterior"))
		  {
		     breadposterior = true;
		  }
		  else if ((args[nargindex].equals("-readstatesbyline"))||(args[nargindex].equals("-readstatebyline")))
		  {
		     breadstatebyline = true;
		  }
                  else if (args[nargindex].equals("-scalebeta"))
	          {
		      bscalebeta = true;
		  }
		  else
		  { 
		     bok = false;
		     break;
		  }
		  nargindex++;
	       }
	    }
	    catch (NumberFormatException ex)
	    {
		bok = false;
	    }

	    if (bok&&(nargindex==args.length-5))
	    {
	       String szmodelfile = args[nargindex++];
	       String szinputdir = args[nargindex++];
	       String szsegmentdir = args[nargindex++];
	       String szconfusionfileprefix = args[nargindex++];
	       String szinclude = args[nargindex];

	       boolean breadsegments = !breadstatebyline&&!breadposterior;
	       if (breadposterior && breadstatebyline)
	       {
		   System.out.println("Invalid to specify both -readposterior and -readstatesbyline output");
	       }
	       else
	       {
		   ChromHMM theHMM = new ChromHMM(szinputdir, szsegmentdir,szinputfilelist,szconfusionfileprefix, 
                                                  szmodelfile, szoutfileID, nbinsize, breadposterior,
						  breadsegments,breadstatebyline,szinclude,bappend, theColor,bprintimage,blowmem, 
                                                  bscaleemissions,bscalebeta);

		  if (blowmem)
		  {
	             theHMM.makeSegmentationConfusionWithLoad();
		  }
		  else
		  {
		     theHMM.makeSegmentationConfusion();
		  }
	       }

	    }
	    else
	    {
		bok = false;
	    }

	    if (!bok)
	    {
		System.out.println("usage: EvalSubset [-append][-b binsize][-f inputfilelist][-i outfileID]"+
                                   "[-lowmem][-many][-noimage][-readposterior|-readstatesbyline][-scalebeta]"+
                                   "  inputmodel inputdir segmentdir outconfusionfileprefix includemarks");
	    }
        }
	else if (szcommand.equalsIgnoreCase("MakeSegmentation"))
	{
	    String szinputfilelist = null;
            boolean bscaleemissions = false;
	    boolean bprintposterior = false;
	    boolean bprintstatebyline = false;
	    boolean bnoprintsegment = false;
	    boolean bscalebeta = false;
	    String szchromlengthfile = null;
	    int nbinsize = ChromHMM.DEFAULT_BINSIZEBASEPAIRS;
	    //int numsplitbins = ChromHMM.DEFAULT_NUMSPLITBINS;
	    String szoutfileID = "";
	    boolean blowmem = false;
	    boolean bgzip = false;
	    int nargindex = 1;
	    boolean bsplit = false;

            try
	    {
	       while (nargindex < args.length-3)
	       {
	          if (args[nargindex].equals("-b"))
		  {
		     nbinsize = Integer.parseInt(args[++nargindex]);
		  }
		  else if (args[nargindex].equals("-f"))
		  {
		     szinputfilelist = args[++nargindex];
		  }
                  else if (args[nargindex].equals("-gzip"))
		  {
		      bgzip = true;
		  }
		  else if (args[nargindex].equals("-i"))
		  {
		     szoutfileID = args[++nargindex];
		  }
		  //else if (args[nargindex].equals("-j"))
	          //{
	       	  //   numsplitbins = Integer.parseInt(args[++nargindex]);
		  //   bsplit = true;
		  //}
		  else if (args[nargindex].equals("-l"))
		  {
		     szchromlengthfile = args[++nargindex];
		  }
		  else if (args[nargindex].equals("-lowmem"))
		  {
		      blowmem = true;
		  }
                  else if (args[nargindex].equals("-many"))
	          {
	             bscaleemissions = true;
		  }
		  else if (args[nargindex].equals("-nobed"))
		  {
		     bnoprintsegment = true;
		  }
		  else if (args[nargindex].equals("-printposterior"))
		  {
		     bprintposterior = true;
		  }
		  else if ((args[nargindex].equals("-printstatesbyline"))||(args[nargindex].equals("-printstatebyline")))
		  {
		     bprintstatebyline = true;
		  }
                  else if (args[nargindex].equals("-scalebeta"))
	          {
		     bscalebeta = true;
		  }
                  else if (args[nargindex].equals("-splitrows"))
	          {
		     bsplit = true;
		  }
		  else
		  { 
		     bok = false;
		     break;
		  }
		  nargindex++;
	       }
	    }
	    catch (NumberFormatException ex)
	    {
		bok = false;
	    }

	    if (bok&&(nargindex==args.length-3))
	    {
	       String szmodelfile = args[nargindex++];
	       String szinputdir = args[nargindex++];
	       String szoutputdir = args[nargindex];
	       File f = new File(szoutputdir);
	       if (!f.exists())
	       {
	          if (!f.mkdirs())
	          {
	             throw new IllegalArgumentException(szoutputdir+" does not exist and could not be created!");
	          }
	       }
		 
	       if (bprintposterior)
	       {
		   File fposterior = new File(szoutputdir+"/POSTERIOR");
	           if (!fposterior.exists())
	           {
	              if (!fposterior.mkdirs())
	              {
	                 throw new IllegalArgumentException(szoutputdir+"POSTERIOR does not exist and could not be created!");
		      }
		   }
	       }//update in v1.20 to allow creating STATEBYLINE directory even if print posterior not asked for		      

	       if (bprintstatebyline)
	       {
	          File fstatebyline = new File(szoutputdir+"/STATEBYLINE");
	          if (!fstatebyline.exists())
	          {
	             if (!fstatebyline.mkdirs())
	             {
	                throw new IllegalArgumentException(szoutputdir+" STATEBYLINE does not exist and could not be created!");
		     }			  
		  }
	       }

	       boolean bprintsegments = !bnoprintsegment;
	       if (bprintsegments||bprintposterior||bprintstatebyline)
	       {

		   ChromHMM theHMM = new ChromHMM(szinputdir, szinputfilelist,szchromlengthfile, szoutputdir, szmodelfile, szoutfileID, nbinsize, bprintposterior,
						  bprintsegments,bprintstatebyline, blowmem,bscaleemissions, bgzip, bsplit, bscalebeta);

		  if (blowmem)
		  {
		     theHMM.makeSegmentationWithLoad();
		  }
		  else
		  {
	             theHMM.makeSegmentation();
		  }
	       }
	       else
	       {
		   System.out.println("No output type was requested!");
	       }
	    }
	    else
	    {
		bok = false;
	    }

	    if (!bok)
	    {
		System.out.println("usage: MakeSegmentation [-b binsize][-f inputfilelist][-gzip][-i outfileID][-l chromosomelengthfile][-lowmem][-many][-nobed]"+
                                   "[-printposterior][-printstatebyline][-scalebeta][-splitrows]"+
                                   "  modelfile inputdir outputdir");
	    }
	}
	else if (szcommand.equalsIgnoreCase("MakeBrowserFiles"))
	{
	    String szcolormapping = null;
	    String szlabelmapping = null;
	    boolean bgzip = false;
	    boolean blowmem = false;
	    boolean bnobrowserheader = false;
	    int nargindex = 1;
	    int numstates = -1;           

	    while (nargindex < args.length-3)
	    {
	       if (args[nargindex].equals("-c"))
	       {
		   szcolormapping = args[++nargindex];
	       }
	       else if (args[nargindex].equals("-gzip"))
	       {
		   bgzip = true;
	       }
	       else if (args[nargindex].equals("-nobrowserheader"))
	       {
		   bnobrowserheader = true;
	       }
	       else if ((args[nargindex].equals("-l"))||(args[nargindex].equals("-m")))
	       {
		   //the -l is for backwards compatibility
		   szlabelmapping = args[++nargindex];
	       }
	       else if (args[nargindex].equals("-lowmem"))
	       {
	          blowmem = true;
	       }
               else if (args[nargindex].equals("-n"))
	       {
	          numstates = Integer.parseInt(args[++nargindex]);
	       }
	       else
	       {
		   bok = false;
		   break;
	       }
	       nargindex++;
	    }

	    if ((bok)&&(nargindex==args.length-3))
	    {
	       String szsegmentfile = args[nargindex++];
	       String szsegmentationname = args[nargindex++];
	       String szoutputfileprefix =args[nargindex];
	       BrowserOutput theBrowserOutput = new BrowserOutput(szsegmentfile,szcolormapping,szlabelmapping,
                                                                  szsegmentationname, szoutputfileprefix,numstates,bgzip, bnobrowserheader);
	       theBrowserOutput.makebrowserdense();
	       if (blowmem)
	       {
	          theBrowserOutput.makebrowserexpandedLowMem();
	       }
	       else
	       {
	          theBrowserOutput.makebrowserexpanded();
	       }
	    }
	    else
	    {
		bok = false;
	    }
	    
	    if (!bok)
	    {
		System.out.println("usage: MakeBrowserFiles [-c colormappingfile][-gzip][-lowmem][-m labelmappingfile][-n numstates][-nobrowserheader] segmentfile segmentationname outputfileprefix");
	    }
	}
	else if (szcommand.equalsIgnoreCase("OverlapEnrichment"))
	{	    
	    String szinput;
	    String szcell = "";
	    String szinputcoorddir;
	    String szinputcoordlist=null;
	    String szlabelmapping = null;
	    int noffsetleft = ChromHMM.DEFAULT_OVERLAPENRICHMENT_NOFFSETLEFT;   //int ChromHMM.DEFAULT_OVERLAPENRICHMENT_NOFFSETLEFT = 0;
	    int noffsetright = ChromHMM.DEFAULT_OVERLAPENRICHMENT_NOFFSETRIGHT;  //int ChromHMM.DEFAULT_OVERLAPENRICHMENT_NOFFSETRIGHT = 1;
	    int nbinsize = ChromHMM.DEFAULT_BINSIZEBASEPAIRS; //
	    boolean bcenter = ChromHMM.DEFAULT_OVERLAPENRICHMENT_BCENTER; //boolean ChromHMM.DEFAULT_OVERLAPENRICHMENT_BCENTER = false;
	    boolean bcountmulti = ChromHMM.DEFAULT_OVERLAPENRICHMENT_BCOUNTMULTI;  //boolean ChromHMM.DEFAULT_OVERLAPENRICHMENT_BCOUNTMULTI = false;
	    boolean busesignal = ChromHMM.DEFAULT_OVERLAPENRICHMENT_BUSESIGNAL;  //boolean ChromHMM.DEFAULT_OVERLAPENRICHMENT_BUSESIGNAL = false;
	    String szcolfields = null; 
	    boolean bbaseres = ChromHMM.DEFAULT_OVERLAPENRICHMENT_BBASERES; //boolean ChromHMM.DEFAULT_OVERLAPENRICHMENT_BBASERES = true;
	    String szoutfile;
	    boolean buniformheat = ChromHMM.DEFAULT_OVERLAPENRICHMENT_BUNIFORMHEAT; //boolean ChromHMM.DEFAULT_OVERLAPENRICHMENT_BUNIFORMHEAT = false;
	    boolean bstringlabels = false;
	    boolean bbrowser = false;


	    String sztitle = "Fold Enrichments";
	    
	    boolean bmax= true;
	    boolean blowmem = false;
	    boolean bprintimage = true;

	    boolean blogpvals = true;

	    int nr=ChromHMM.DEFAULTCOLOR_R;
	    int ng=ChromHMM.DEFAULTCOLOR_G;
	    int nb=ChromHMM.DEFAULTCOLOR_B;

	    int nargindex = 1;

            try
	    {
	       while (nargindex < args.length-3)
	       {
		  if (args[nargindex].equals("-a"))
		  {
		     szcell = args[++nargindex];
		  }
		  else if (args[nargindex].equals("-b"))
		  {
		     nbinsize = Integer.parseInt(args[++nargindex]);
		  }
		  else if (args[nargindex].equals("-color"))
		  {
		     String szcolor = args[++nargindex];
		     StringTokenizer stcolor = new StringTokenizer(szcolor,",");
		     if (stcolor.countTokens()==3)
		     {
			nr = Integer.parseInt(stcolor.nextToken().trim());
			ng = Integer.parseInt(stcolor.nextToken().trim());
			nb = Integer.parseInt(stcolor.nextToken().trim());
		     }
		     else
		     {
			bok = false;
		     }
		  }
		  else if (args[nargindex].equals("-colfields"))
		  {
		     szcolfields = args[++nargindex];
		  }
		  else if (args[nargindex].equals("-binres"))
		  {
		     bbaseres = false;
		  }
		  else if (args[nargindex].equals("-browser"))
		  {
		      bbrowser = true;
		  }
		  else if (args[nargindex].equals("-f"))
		  {
		     szinputcoordlist = args[++nargindex];
		  }
		  else if (args[nargindex].equals("-nopvals"))
		  {
		      blogpvals = false;
	          }
		  else if (args[nargindex].equals("-signal"))
		  {
		     busesignal = true;
		  }
		  else if (args[nargindex].equals("-uniformscale"))
		  {
		     buniformheat = true; //scales heatmap columns individually
		  }
		  else if (args[nargindex].equals("-s"))
		  {
		     noffsetleft = Integer.parseInt(args[++nargindex]);
		  }
		  else if (args[nargindex].equals("-center"))
		  {
		     bcenter = true;
		  }
		  else if (args[nargindex].equals("-labels"))
		  {
		      bstringlabels = true;
		  }
                  else if (args[nargindex].equals("-m"))
		  {
		      szlabelmapping = args[++nargindex];
		  }
		  else if (args[nargindex].equals("-multicount"))
		  {
		     bcountmulti = true;
		  }
                  else if (args[nargindex].equals("-noimage"))
		  {
		     bprintimage = false;
	          }
		  else if (args[nargindex].equals("-posterior"))
		  {
		     bmax = false;
		  }
		  else if (args[nargindex].equals("-e"))
		  {
		     noffsetright = Integer.parseInt(args[++nargindex]);
		  }
		  else if (args[nargindex].equals("-t"))
		  {
		     sztitle = args[++nargindex];
		  }
		  else if (args[nargindex].equals("-lowmem"))
		  {
		      blowmem = true;
		  }
		  else
	          { 
		     bok = false;
		     break;
		  }
		  nargindex++;
	       }
	    }
	    catch (NumberFormatException ex)
	    {
		bok = false;
	    }

	    if ((nargindex == args.length-3)&&(bok))
	    {
	       szinput = args[nargindex++];
	       szinputcoorddir = args[nargindex++];
	       szoutfile = args[nargindex];
	       boolean bunique = !bcountmulti;
	       boolean bcolscaleheat = !buniformheat;
	       Color theColor = new Color(nr, ng, nb);
	       System.out.println("Computing Enrichments...");
	       if (bmax)
	       {
		   if (blowmem)
		   {
		      StateAnalysis.enrichmentMaxLowMem(szinput, szinputcoorddir,szinputcoordlist,noffsetleft,noffsetright,nbinsize,  
							bcenter, bunique,  busesignal,szcolfields,bbaseres, szoutfile,bcolscaleheat,theColor,sztitle, szlabelmapping, bprintimage,bstringlabels,bbrowser,blogpvals);
		   }
		   else
		   {
		      StateAnalysis.enrichmentMax(szinput, szinputcoorddir,szinputcoordlist,noffsetleft,noffsetright,nbinsize,  
						  bcenter, bunique,  busesignal,szcolfields,bbaseres, szoutfile,bcolscaleheat,theColor,sztitle, szlabelmapping,bprintimage,bstringlabels,bbrowser,blogpvals);
		   }
	       }
	       else
	       {
		   StateAnalysis.enrichmentPosterior(szinput, szcell,szinputcoorddir,szinputcoordlist,noffsetleft,noffsetright,nbinsize,
						     bcenter, bunique, busesignal,szcolfields,bbaseres,szoutfile,bcolscaleheat,theColor,sztitle, szlabelmapping,bprintimage,blogpvals);
	       }
	    }
	    else
	    {
		bok = false;
	    }

	    if (!bok)
	    {
		System.out.println("usage OverlapEnrichment [-a cell][-b binsize][-binres][-browser][-color r,g,b][-center][-colfields chromosome,start,end[,signal]]"+
                                   "[-e offsetend][-f coordlistfile][-labels][-lowmem][-m labelmappingfile]"+
                                   "[-multicount][-noimage][-nopvals][-posterior][-s offsetstart][-signal][-t title][-uniformscale]"+
                                   " inputsegment inputcoorddir outfileprefix"); 
	    }
	}
	else if ((szcommand.equalsIgnoreCase("NeighborhoodEnrichment"))||(szcommand.equalsIgnoreCase("NeighborhoodSignal")))
	{
	    int nbinsize = ChromHMM.DEFAULT_BINSIZEBASEPAIRS;
	    String szinput;

	    int numleft = ChromHMM.DEFAULT_NEIGHBORHOOD_NUMLEFT; 
	    int numright = ChromHMM.DEFAULT_NEIGHBORHOOD_NUMRIGHT; 
	    int nspacing = ChromHMM.DEFAULT_BINSIZEBASEPAIRS;
	    String szlabelmapping = null;
	    boolean busestrand =ChromHMM.DEFAULT_NEIGHBORHOOD_BUSESTRAND; 
	    boolean busesignal = ChromHMM.DEFAULT_NEIGHBORHOOD_BUSESIGNAL; 
	    boolean bstringlabels = false;
	    boolean bbrowser = false;
	    String szcolfields = null;
	    int noffsetanchor = ChromHMM.DEFAULT_NEIGHBORHOOD_NOFFSETANCHOR; 
	    String szoutfile; 
	    String sztitle = "Fold Enrichments";
	    String szcell = ""; 
	    String szanchorpositions; 

	    int nr=ChromHMM.DEFAULTCOLOR_R;
	    int ng=ChromHMM.DEFAULTCOLOR_G;
	    int nb=ChromHMM.DEFAULTCOLOR_B;

	    boolean bmax = true;
	    boolean blowmem = false;
	    boolean bprintimage = true;

	    int nargindex = 1;

	    boolean bspacing = false;

            try
	    {
	       while (nargindex < args.length-3)
	       {
		  if (args[nargindex].equals("-a"))
		  {
		     szcell = args[++nargindex];
		  }
		  else if (args[nargindex].equals("-b"))
		  {
		     nbinsize = Integer.parseInt(args[++nargindex]);
		  }
		  else if (args[nargindex].equals("-browser"))
		  {
		      bbrowser = true;
		  }
		  else if (args[nargindex].equals("-color"))
		  {
		     String szcolor = args[++nargindex];
		     StringTokenizer stcolor = new StringTokenizer(szcolor,",");
		     if (stcolor.countTokens()==3)
		     {
			nr = Integer.parseInt(stcolor.nextToken().trim());
			ng = Integer.parseInt(stcolor.nextToken().trim());
			nb = Integer.parseInt(stcolor.nextToken().trim());
		     }
		     else
		     {
			bok = false;
		     }
		  }
		  else if (args[nargindex].equals("-colfields"))
		  {
		     szcolfields = args[++nargindex];
		  }
		  else if (args[nargindex].equals("-signal"))
		  {
		     busesignal = true;
		  }
		  else if (args[nargindex].equals("-l"))
		  {
		     numleft = Integer.parseInt(args[++nargindex]);
		  }
		  else if (args[nargindex].equals("-labels"))
		  {
		      bstringlabels = true;
		  }
		  else if (args[nargindex].equals("-lowmem"))
		  {
		      blowmem = true;
		  }
		  else if (args[nargindex].equals("-m"))
	          {
	             szlabelmapping = args[++nargindex];
		  }
                  else if (args[nargindex].equals("-noimage"))
		  {
		      bprintimage = false;
		  }
		  else if (args[nargindex].equals("-o"))
		  {
		     noffsetanchor = Integer.parseInt(args[++nargindex]);
		  }
		  else if (args[nargindex].equals("-posterior"))
		  {
		     bmax = false;
		  }
		  else if (args[nargindex].equals("-r"))
		  {
		     numright = Integer.parseInt(args[++nargindex]);
		  }
		  else if (args[nargindex].equals("-s"))
		  {
		     nspacing = Integer.parseInt(args[++nargindex]);
		     bspacing = true;
		  }
		  else if (args[nargindex].equals("-nostrand"))
		  {
		     busestrand = false;
		  }
		  else if (args[nargindex].equals("-t"))
		  {
		     sztitle = args[++nargindex];
		  }
		  else
		  {
		     bok =false;
		     break;
		  }
		  nargindex++;
	       }
	    }
	    catch (NumberFormatException ex)
	    {
		bok = false;
	    }


	    if (!bspacing)
	    {
		nspacing = nbinsize;
	    }

	    if ((nargindex == args.length-3)&&(bok))
	    {
	       Color theColor = new Color(nr, ng, nb);
	       szinput = args[nargindex++];
	       szanchorpositions = args[nargindex++];
	       szoutfile = args[nargindex];
	       System.out.println("Computing Enrichments...");
	       if (szcommand.equalsIgnoreCase("NeighborhoodSignal"))
	       {
		   //this is an undocumented feature to compute signal enrichment for marks around a position
	           StateAnalysis.neighborhoodSignal(szinput,szcell,szanchorpositions,nbinsize,numleft,numright,
						    nspacing,busestrand,busesignal,szcolfields,noffsetanchor,szoutfile,theColor,
                                                    sztitle,szlabelmapping, bprintimage);
	       }
               else if (bmax)
	       {
		   if (blowmem)
		   {
	              StateAnalysis.neighborhoodMaxLowMem(szinput,szanchorpositions,nbinsize,numleft,numright,
							  nspacing,busestrand,busesignal,szcolfields,noffsetanchor,szoutfile,theColor,
                                                          sztitle,szlabelmapping, bprintimage, bstringlabels, bbrowser);
		   }
		   else
		   {
	              StateAnalysis.neighborhoodMax(szinput,szanchorpositions,nbinsize,numleft,numright,
						    nspacing,busestrand,busesignal,szcolfields,noffsetanchor,szoutfile,theColor,sztitle,szlabelmapping, 
                                                    bprintimage, bstringlabels,bbrowser);
		   }
	       }
	       else
	       {
	           StateAnalysis.neighborhoodPosterior(szinput,szcell,szanchorpositions,nbinsize,numleft,numright,
						       nspacing,busestrand,busesignal,szcolfields,noffsetanchor,szoutfile,theColor,sztitle,szlabelmapping,bprintimage);
	       }
	    }
	    else
	    {
	       bok = false;
	    }

	    if (!bok)
	    {
		System.out.println("usage NeighborhoodEnrichment [-a cell][-b binsize][-browser][-color r,g,b][-colfields chromosome,position[,optionalcol1|,optionalcol1,optionalcol2]"+
                                  "[-l numleftintervals][-labels][-lowmem][-m labelmappingfile][-noimage][-nostrand]"+
                                  "[-o anchoroffset][-posterior][-r numrightintervals]"+
                                  "[-s spacing][-signal][-t title] inputsegment anchorpositions outfileprefix");
	    }
	}
	else if (szcommand.equalsIgnoreCase("LearnModel"))
        {
	    String path = ChromHMM.class.getProtectionDomain().getCodeSource().getLocation().getPath();
	    String decodedPath = URLDecoder.decode(path, "UTF-8");
	    String szprefixpath = decodedPath.substring(0, decodedPath.lastIndexOf("/") + 1); 
	    //String szfullprefixpathchromsizedir = szprefixpath+"/"+CHROMSIZESDIR;
	    String szfullprefixpathcoorddir = szprefixpath+"/"+COORDDIR;
	    String szfullprefixpathanchorfiledir = szprefixpath+"/"+ANCHORFILEDIR;

	    //System.out.println("prefix path is "+szprefixpath);


	    int nseed = 999;
	    int nmaxiterations= 200;//number of passes through all the data
	    int nzerotransitionpower =8; //Sets the transition probability to 0 if below this constant value
	    String szInitFile =null;
	    int ninitmethod = ChromHMM.INITMETHOD_INFORMATION;
	    double dconvergediff = 0.001;//-1;
	    double dinformationsmooth= 0.02;
	    double dloadsmoothemission = 0.02;
	    double dloadsmoothtransition = 0.5;
	    boolean bautoopen = true;
	    boolean bprintposterior = false;
	    boolean bprintstatebyline = false;
	    boolean bnoprintsegment = false;
	    boolean bprintbrowser = true;
	    boolean bprintenrich = true;
	    boolean bnobrowserheader = false;
	    boolean bprintimage = true;
	    boolean blowmem = false;
	    boolean bscaleemissions = false;
	    boolean bpseudo = false;
	    boolean bnopseudolast = false;
	    boolean bgzip = false;
	    boolean bsplit = false;
	    int numsplitbins = ChromHMM.DEFAULT_NUMSPLITBINS;

	    String szinputfilelist = null;
	    String szchromlengthfile = null;
	    int nbinsize = ChromHMM.DEFAULT_BINSIZEBASEPAIRS;
	    String szoutfileID = "";
	    int nstateorder =ChromHMM.STATEORDER_EMISSION;
	    boolean bnoordercols = false;
	    boolean bnoorderrows = false;
	    int nargindex = 1;
	    int nmaxseconds = -1;
	    int nmaxprocessors = 0;
	    int numincludeseq = 0;
	    boolean bnormalEM = false;
	    boolean bscalebeta = false;

	    int nr=ChromHMM.DEFAULTCOLOR_R;
	    int ng=ChromHMM.DEFAULTCOLOR_G;
	    int nb=ChromHMM.DEFAULTCOLOR_B;

            try
	    {
	       while (nargindex < args.length-4)
	       {
		  if (args[nargindex].equals("-b"))
		  {
		     nbinsize = Integer.parseInt(args[++nargindex]);
		  }
		  else if (args[nargindex].equals("-color"))
		  {
		     String szcolor = args[++nargindex];
		     StringTokenizer stcolor = new StringTokenizer(szcolor,",");
		     if (stcolor.countTokens()==3)
		     {
		        nr = Integer.parseInt(stcolor.nextToken().trim());
		        ng = Integer.parseInt(stcolor.nextToken().trim());
		        nb = Integer.parseInt(stcolor.nextToken().trim());
		     }
		     else
		     {
		        bok = false;
		     }
		  }
		  else if (args[nargindex].equals("-d"))
		  {
		     dconvergediff = Double.parseDouble(args[++nargindex]);
		  }
		  else if (args[nargindex].equals("-e"))
		  {
		     dloadsmoothemission = Double.parseDouble(args[++nargindex]);
		  }
		  else if (args[nargindex].equals("-f"))
		  {
		     szinputfilelist = args[++nargindex];
		  }
                  else if (args[nargindex].equals("-gzip"))
		  {
	       	     bgzip = true;
		  }
		  else if (args[nargindex].equals("-h"))
		  {
		      dinformationsmooth = Double.parseDouble(args[++nargindex]);
		  }
		  else if (args[nargindex].equals("-i"))
		  {
		     szoutfileID = args[++nargindex];
		  }
		  else if (args[nargindex].equals("-init"))
		  {
                     String sztoken = args[++nargindex];

		     if ((sztoken.equals("information"))||(sztoken.equals("0")))
		     {
			 ninitmethod = ChromHMM.INITMETHOD_INFORMATION;
		     }
		     else if ((sztoken.equals("random"))||(sztoken.equals("1")))
		     {
			 ninitmethod = ChromHMM.INITMETHOD_RANDOM;
		     } 
                     else if ((sztoken.equals("load"))||(sztoken.equals("2")))
		     {
			 ninitmethod = ChromHMM.INITMETHOD_LOAD;
		     }
		     else
		     {
			 bok = false;
			 break;
		     }
		  }
		  //else if (args[nargindex].equals("-j"))
	          //{
                  //   numsplitbins = Integer.parseInt(args[++nargindex]);
		  //   bsplit = true;
		  // }
		  else if (args[nargindex].equals("-l"))
		  {
		      szchromlengthfile = args[++nargindex];
		  }
		  else if (args[nargindex].equals("-m"))
		  {
		     szInitFile = args[++nargindex];
		  }
                  else if (args[nargindex].equals("-many"))
	          {
		      bscaleemissions = true;
		  }
		  else if (args[nargindex].equals("-n"))
		  {
		      numincludeseq = Integer.parseInt(args[++nargindex]);
		  }
                  else if (args[nargindex].equals("-noautoopen"))
		  {
		     bautoopen = false;
		  }
		  else if (args[nargindex].equals("-nobed"))
		  {
		     bnoprintsegment = true;
		  }
		  else if (args[nargindex].equals("-nobrowser"))
		  {
		      bprintbrowser = false;
		  }
		  else if (args[nargindex].equals("-nobrowserheader"))
		  {
		      bnobrowserheader = true;
		  }
		  else if (args[nargindex].equals("-noenrich"))
		  {
		      bprintenrich = false;
		  }
                  else if (args[nargindex].equals("-noimage"))
		  {
		      bprintimage = false;
		  }
		  else if (args[nargindex].equals("-p"))
		  {
		      bnormalEM = true;
		      nmaxprocessors = Integer.parseInt(args[++nargindex]);
		  }
                  else if (args[nargindex].equals("-pseudo"))
		  {
		     bpseudo = true;
		  }
		  else if (args[nargindex].equals("-nopseudolast"))
		  {
		      bnopseudolast = true;
		  }
		  else if (args[nargindex].equals("-stateordering"))
		  {
                     String sztoken = args[++nargindex];
                     if (sztoken.equals("emission"))
		     {
			 nstateorder = ChromHMM.STATEORDER_EMISSION;
		     } 
                     else if (sztoken.equals("transition"))
		     {
			 nstateorder = ChromHMM.STATEORDER_TRANSITION;
		     }
		     else
		     {
			 bok = false;
			 break;
		     }
		  }
		  else if (args[nargindex].equals("-printposterior"))
		  {
		     bprintposterior = true;
		  }
                  else if (args[nargindex].equals("-lowmem"))
	          {
		      blowmem = true;
		  }
		  else if (args[nargindex].equals("-r"))
		  {
		     nmaxiterations = Integer.parseInt(args[++nargindex]);
		  }	
		  else if (args[nargindex].equals("-s"))
		  {
		     nseed = Integer.parseInt(args[++nargindex]);
		  }
                  else if (args[nargindex].equals("-scalebeta"))
		  {
		      bscalebeta = true;
		  }
		  else if (args[nargindex].equals("-t"))
		  {
		     dloadsmoothtransition = Double.parseDouble(args[++nargindex]);
		  }
		  else if (args[nargindex].equals("-holdcolumnorder"))
		  {
		     bnoordercols = true;
		  }
                  else if (args[nargindex].equals("-holdroworder"))
	          {
		     bnoorderrows = true;
		     nstateorder = ChromHMM.STATEORDER_FIXED;
		  }
		  else if ((args[nargindex].equals("-printstatebyline"))||(args[nargindex].equals("-printstatesbyline")))
		  {
		      bprintstatebyline = true;
		  }
                  else if (args[nargindex].equals("-splitrows"))
	          {
		      bsplit = true;
		  }
		  //else if (args[nargindex].equals("-u"))
		  //{
		  //    szfullprefixpathchromsizedir = args[++nargindex];
		  //}
		  else if (args[nargindex].equals("-u"))
		  {
		      szfullprefixpathcoorddir = args[++nargindex];
		  }
                  else if (args[nargindex].equals("-v"))
		  {
	       	     szfullprefixpathanchorfiledir = args[++nargindex];
		  }
		  else if (args[nargindex].equals("-x"))
		  {
		      nmaxseconds = Integer.parseInt(args[++nargindex]);
		  }
		  else if (args[nargindex].equals("-z"))
		  {
		     nzerotransitionpower = Integer.parseInt(args[++nargindex]);
		  }
		  else
		  {
		     bok = false;
		  }
		  nargindex++;
	       }
	    }
	    catch (NumberFormatException ex)
	    {
		bok = false;
	    }


	    if (bok&&(nargindex==args.length-4))
	    {
	       Color theColor = new Color(nr, ng, nb);
	       String szinputdir = args[nargindex++];
	       String szoutputdir = args[nargindex++];
	       File f = new File(szoutputdir);
	       if (!f.exists())
	       { 
	          if (!f.mkdirs())
	          {
	             throw new IllegalArgumentException(szoutputdir+" does not exist and could not be created!");
	          }
	       }

	       if (bprintposterior)
	       {
	          File fposterior = new File(szoutputdir+"/POSTERIOR");
	          if (!fposterior.exists())
	          {
	             if (!fposterior.mkdirs())
	             {
	                throw new IllegalArgumentException(szoutputdir+"/POSTERIOR does not exist and could not be created!");
		     }
		  }
	       }   

	       if (bprintstatebyline)
	       {
		   File fstatebyline = new File(szoutputdir+"/STATEBYLINE");
	           if (!fstatebyline.exists())
	           {
	              if (!fstatebyline.mkdirs())
	              {
	                 throw new IllegalArgumentException(szoutputdir+"/STATEBYLINE does not exist and could not be created!");		           
		      }
		   }		  
	       }

	       int numstates = 0;
               try
	       {
	          numstates = Integer.parseInt(args[nargindex++]);
	       }
	       catch (NumberFormatException ex)
	       {
		   bok = false;
	       }

	       if (bok)
	       {
		  String szassembly= args[nargindex++];

	          boolean bprintsegments = !bnoprintsegment; 
	          boolean bordercols = !bnoordercols;
		  boolean borderrows = !bnoorderrows;
		  if ((szoutfileID.equals(""))&&(ninitmethod == ChromHMM.INITMETHOD_RANDOM))
		  {
		      szoutfileID = ""+nseed;
		  }

		  if (szchromlengthfile == null)
		  {
		      //updated v.1.16
		      File flength = new File(szprefixpath+"/"+CHROMSIZESDIR+"/"+szassembly+".txt");
		      //File flength = new File(szfullprefixpathchromsizedir+"/"+szassembly+".txt");
		      if (flength.exists())
		      {
			  //szchromlengthfile = szfullprefixpathchromsizedir+"/"+szassembly+".txt";
			  szchromlengthfile = szprefixpath+"/"+CHROMSIZESDIR+"/"+szassembly+".txt";
		      }
		  }
 	          ChromHMM theHMM = new ChromHMM(szinputdir, szoutputdir, szinputfilelist,szchromlengthfile,numstates, nseed,ninitmethod,
						 szInitFile,dloadsmoothemission,dloadsmoothtransition,dinformationsmooth,
					         nmaxiterations,dconvergediff,nmaxseconds, bprintposterior,bprintsegments,bprintstatebyline,
						 nbinsize,szoutfileID,nstateorder,bordercols,nzerotransitionpower,theColor,bnormalEM, nmaxprocessors, 
                                                 blowmem,numincludeseq,bprintimage,bscaleemissions, bpseudo,bnopseudolast, bgzip, bsplit, borderrows, bscalebeta);
	          theHMM.buildModel();


		  String szunderscoreoutfileID;
		  if (szoutfileID.equals(""))
		  {
		      szunderscoreoutfileID = "";
		  }
		  else
		  {
                      szunderscoreoutfileID = "_"+szoutfileID;
		  }

		  String szwebpage = szoutputdir+"/webpage_"+numstates+szunderscoreoutfileID+".html";
		  PrintWriter pwweb = new PrintWriter(new FileWriter(szwebpage));


		  pwweb.println("<center><h1>ChromHMM Report</h1></center>");
		  pwweb.println("Input Directory: "+szinputdir+"<br>");
                  pwweb.println("Output Directory: "+szoutputdir+"<br>");
                  pwweb.println("Number of States: "+numstates+"<br>");
                  pwweb.println("Assembly: "+szassembly+"<br>");
		  pwweb.print("Full ChromHMM command: ");
		  for (int na = 0; na < args.length-1; na++)
		  {
		      pwweb.print(args[na]+" ");
		  }
		  pwweb.println(args[args.length-1]+"");



		  pwweb.println("<h1>Model Parameters</h1>");
		  if (bprintimage)
		  {
		     pwweb.println("<img src=\"emissions_"+numstates+szunderscoreoutfileID+".png\"><br>");
		     pwweb.println("<li><a href=\"emissions_"+numstates+szunderscoreoutfileID+".svg\">Emission Parameter SVG File</a><br>");
		  }
		  pwweb.println("<li><a href=\"emissions_"+numstates+szunderscoreoutfileID+".txt\">Emission Parameter Tab-Delimited Text File</a><br>");
		  if (bprintimage)
		  {
		     pwweb.println("<img src=\"transitions_"+numstates+szunderscoreoutfileID+".png\"><br>");
	             pwweb.println("<li><a href=\"transitions_"+numstates+szunderscoreoutfileID+".svg\">Transition Parameter SVG File</a><br>");
		  }
		  pwweb.println("<li><a href=\"transitions_"+numstates+szunderscoreoutfileID+".txt\">Transition Parameter Tab-Delimited Text File</a><br><br>");
		  pwweb.println("<li><a href=\"model_"+numstates+szunderscoreoutfileID+".txt\">All Model Parameters Tab-Delimited Text File</a> <br>");
		  pwweb.println("<h1>Genome Segmentation Files</h1>");

	          if (bprintsegments||bprintposterior||bprintstatebyline) 
	          {
		     if (blowmem)
		     {
			 theHMM.makeSegmentationWithLoad();
		     }
		     else
		     {
	                theHMM.makeSegmentation();
		     }

		     if (bprintsegments)
		     {
			Iterator hsiterator = theHMM.hsprefix.iterator();
			String[] prefixes = new String[theHMM.hsprefix.size()];
			int nhsindex = 0;
			while (hsiterator.hasNext())
			{
			    prefixes[nhsindex] = (String) hsiterator.next();
			    nhsindex++;
			}
			Arrays.sort(prefixes);

			for (int nfile = 0; nfile < prefixes.length; nfile++)
			{
			   String szsegmentfile;
			   if (bgzip)
			   {
                              szsegmentfile = prefixes[nfile]+ChromHMM.SZSEGMENTEXTENSION+".gz";
			   }
			   else
			   {
                              szsegmentfile = prefixes[nfile]+ChromHMM.SZSEGMENTEXTENSION;
			   }
			   pwweb.println("<li><a href=\""+szsegmentfile+"\">"+prefixes[nfile]+" Segmentation File (Four Column Bed File)</a><br>");
			}

		        if (bprintstatebyline)
			{
			    pwweb.println("<li><a href=\"STATEBYLINE\"> Directory of Maximum States Assignments Line By Line</a><br>");
			}

			if (bprintposterior)
			{
			    pwweb.println("<li><a href=\"POSTERIOR\"> Directory of Posterior Files</a><br>");
			}


			if (bprintbrowser)
			{
			    pwweb.println("<br>");
			    pwweb.println("Custom Tracks for loading into the <a href=\"http://genome.ucsc.edu\">UCSC Genome Browser</a>:<br>");

			    for (int nprefixindex = 0; nprefixindex < prefixes.length; nprefixindex++)
			    {
			      String szprefix = prefixes[nprefixindex];

			      String szsegmentfile;
			      if (bgzip)
			      {
				  szsegmentfile = szoutputdir+"/"+szprefix+ChromHMM.SZSEGMENTEXTENSION+".gz";
			      }
			      else
			      {
			         szsegmentfile = szoutputdir+"/"+szprefix+ChromHMM.SZSEGMENTEXTENSION;
			      }


			      BrowserOutput theBrowserOutput = new BrowserOutput(szsegmentfile,null,null,
										 szprefix,szoutputdir+"/"+szprefix,numstates,bgzip,bnobrowserheader);

			      theBrowserOutput.makebrowserdense();
			      if (blowmem)
			      {
			         theBrowserOutput.makebrowserexpandedLowMem();
			      }
			      else
                              {
			         theBrowserOutput.makebrowserexpanded();
			      }

			      if (bgzip)
			      {
			         pwweb.println("<li><a href="+szprefix+ChromHMM.SZBROWSERDENSEEXTENSION+".bed.gz>"+szprefix+" Browser Custom Track Dense File</a> <br>");
			         pwweb.println("<li><a href="+szprefix+ChromHMM.SZBROWSEREXPANDEDEXTENSION+".bed.gz>"+szprefix+" Browser Custom Track Expanded File</a><br>");
			      }
			      else
			      {
			         pwweb.println("<li><a href="+szprefix+ChromHMM.SZBROWSERDENSEEXTENSION+".bed>"+szprefix+" Browser Custom Track Dense File</a> <br>");
			         pwweb.println("<li><a href="+szprefix+ChromHMM.SZBROWSEREXPANDEDEXTENSION+".bed>"+szprefix+" Browser Custom Track Expanded File</a><br>");
			      }
			   }
			}


			if (bprintenrich)
		        {
                           pwweb.println("<h1>State Enrichments</h1>");

			    for (int nprefixindex = 0; nprefixindex < prefixes.length; nprefixindex++)
			    {
			       String szprefix = prefixes[nprefixindex];
			       pwweb.println("<h2>"+szprefix+" Enrichments</h2>");
			       String szsegmentfile;
			       if (bgzip)
			       {
			          szsegmentfile = szoutputdir+"/"+szprefix+ChromHMM.SZSEGMENTEXTENSION+".gz";
			       }
			       else
			       {
				  szsegmentfile = szoutputdir+"/"+szprefix+ChromHMM.SZSEGMENTEXTENSION;
			       }
			       //File fcoordassembly = new File(szprefixpath+"/"+COORDDIR+"/"+szassembly);
			       File fcoordassembly = new File(szfullprefixpathcoorddir+"/"+szassembly);
		               if (fcoordassembly.exists())
			       {
				  if (blowmem)
			          {
			             StateAnalysis.enrichmentMaxLowMem(szsegmentfile,szfullprefixpathcoorddir+"/"+szassembly,null,//szinputcoordlist,
                                                          ChromHMM.DEFAULT_OVERLAPENRICHMENT_NOFFSETLEFT,
                                                          ChromHMM.DEFAULT_OVERLAPENRICHMENT_NOFFSETRIGHT,nbinsize,  
							  ChromHMM.DEFAULT_OVERLAPENRICHMENT_BCENTER, !ChromHMM.DEFAULT_OVERLAPENRICHMENT_BCOUNTMULTI, 
                                                          ChromHMM.DEFAULT_OVERLAPENRICHMENT_BUSESIGNAL,null,//szcolfields,
                                                          ChromHMM.DEFAULT_OVERLAPENRICHMENT_BBASERES, szoutputdir+"/"+szprefix+ChromHMM.SZOVERLAPEXTENSION,
								       !ChromHMM.DEFAULT_OVERLAPENRICHMENT_BUNIFORMHEAT,theColor,"Fold Enrichment "+szprefix,null,bprintimage, false,false,false);
				  }
				  else
				  {
				      //update v.1.16
				      //StateAnalysis.enrichmentMax(szsegmentfile,szprefixpath+"/"+COORDDIR+"/"+szassembly,null,//szinputcoordlist,
			             StateAnalysis.enrichmentMax(szsegmentfile,szfullprefixpathcoorddir+"/"+szassembly,null,//szinputcoordlist,
                                                          ChromHMM.DEFAULT_OVERLAPENRICHMENT_NOFFSETLEFT,
                                                          ChromHMM.DEFAULT_OVERLAPENRICHMENT_NOFFSETRIGHT,nbinsize,  
							  ChromHMM.DEFAULT_OVERLAPENRICHMENT_BCENTER, !ChromHMM.DEFAULT_OVERLAPENRICHMENT_BCOUNTMULTI, 
                                                          ChromHMM.DEFAULT_OVERLAPENRICHMENT_BUSESIGNAL,null,//szcolfields,
                                                          ChromHMM.DEFAULT_OVERLAPENRICHMENT_BBASERES, szoutputdir+"/"+szprefix+ChromHMM.SZOVERLAPEXTENSION,
								 !ChromHMM.DEFAULT_OVERLAPENRICHMENT_BUNIFORMHEAT,theColor,"Fold Enrichment "+szprefix,null,bprintimage, false,false,false);
				  }
				  String szoverlapoutfile = szprefix+ChromHMM.SZOVERLAPEXTENSION+".txt";
				  if (bprintimage)
				  {
				     pwweb.println("<img src=\""+szprefix+ChromHMM.SZOVERLAPEXTENSION+".png\"> <br>");
				     pwweb.println("<li><a href=\""+szprefix+ChromHMM.SZOVERLAPEXTENSION+".svg"+"\">"+szprefix+" Overlap Enrichment SVG File"+"</a><br>");
				  }
				  pwweb.println("<li><a href=\""+szoverlapoutfile+"\">"+szprefix+" Overlap Enrichment Tab-Delimited Text File"+"</a><br>");

			       }
			       else
			       {
				   System.out.println("Warning: No coordinate directory found for assembly "+szassembly+" in "+szfullprefixpathcoorddir);//szprefixpath+"/"+COORDDIR);
			       }

                               //szfullprefixpathanchorfiledir
			       //File fanchorassembly = new File(szprefixpath+"/"+ANCHORFILEDIR+"/"+szassembly);
			       File fanchorassembly = new File(szfullprefixpathanchorfiledir+"/"+szassembly);
			       if (fanchorassembly.exists())
			       {
                                  String[] dir = fanchorassembly.list();
			          for (int nfile = 0; nfile < dir.length; nfile++)
			          {
				     File fhidden = new File(dir[nfile]);
				     if (fhidden.isHidden())
				     {
					 //added in v1.11
				        continue;
				     }

				     int nperiodindex = dir[nfile].indexOf(".");
				     String szanchorname;
				     if (nperiodindex == -1)
				     {
				        szanchorname = dir[nfile];
				     }
				     else
				     {
				        szanchorname = dir[nfile].substring(0,nperiodindex);
				     }

				     if (blowmem)
				     {
			                StateAnalysis.neighborhoodMaxLowMem(szsegmentfile,
                                                                 szfullprefixpathanchorfiledir+"/"+szassembly+"/"+dir[nfile],nbinsize,ChromHMM.DEFAULT_NEIGHBORHOOD_NUMLEFT,
							     ChromHMM.DEFAULT_NEIGHBORHOOD_NUMRIGHT, nbinsize,//nspacing
                                           ChromHMM.DEFAULT_NEIGHBORHOOD_BUSESTRAND,ChromHMM.DEFAULT_NEIGHBORHOOD_BUSESIGNAL,null,//szcolfields,
                                           ChromHMM.DEFAULT_NEIGHBORHOOD_NOFFSETANCHOR,szoutputdir+"/"+szprefix+"_"+szanchorname+"_neighborhood",
									    theColor,"Fold Enrichment "+szprefix+" "+szanchorname,null,bprintimage,false,false);
				     }
				     else
				     {
			                StateAnalysis.neighborhoodMax(szsegmentfile,
                                                                 szfullprefixpathanchorfiledir+"/"+szassembly+"/"+dir[nfile],nbinsize,ChromHMM.DEFAULT_NEIGHBORHOOD_NUMLEFT,
							     ChromHMM.DEFAULT_NEIGHBORHOOD_NUMRIGHT, nbinsize,//nspacing
                                           ChromHMM.DEFAULT_NEIGHBORHOOD_BUSESTRAND,ChromHMM.DEFAULT_NEIGHBORHOOD_BUSESIGNAL,null,//szcolfields,
                                           ChromHMM.DEFAULT_NEIGHBORHOOD_NOFFSETANCHOR,szoutputdir+"/"+szprefix+"_"+szanchorname+"_neighborhood",
								      theColor,"Fold Enrichment "+szprefix+" "+szanchorname,null,bprintimage,false,false);
				     }
				     String szneighborhoodoutfileprefix = szprefix+"_"+szanchorname+ChromHMM.SZNEIGHBORHOODEXTENSION;
				   
				     if (bprintimage)
				     {
				        pwweb.println("<img src=\""+szneighborhoodoutfileprefix+".png\"> <br>");
				        pwweb.println("<li><a href=\""+szneighborhoodoutfileprefix+".svg"+"\">"+szneighborhoodoutfileprefix+" Enrichment SVG File</a><br>");
				     }
				     pwweb.println("<li><a href=\""+szneighborhoodoutfileprefix+".txt"+"\">"+szneighborhoodoutfileprefix+" Enrichment Tab-Delimited Text File</a><br>");
				  }			       
			       }
			       else
			       {
			          System.out.println("Warning: No coordinate directory found for assembly "+szassembly+" in "+szfullprefixpathanchorfiledir);
			       }
			   }
			}
		     }
		  }
	          pwweb.close();
                 
		  if (bautoopen)
		  {
                     try
		     {
		        java.awt.Desktop.getDesktop().browse((new File(szwebpage)).toURI());
		     }
		     catch (Exception ex)
		     {
		        System.out.println("Warning could not automatically open in a browser "+szwebpage);
		     }
		  }		 
	       }
	    }
	    else
	    {
		bok = false;
	    }

	    if (!bok)
	    {
		System.out.println("usage: LearnModel [-b binsize][-color r,g,b][-d convergedelta][-e loadsmoothemission][-f inputfilelist][-gzip][-h informationsmooth]"+
                                     "[-holdcolumnorder][-holdroworder][-i outfileID][-init information|random|load][-l chromosomelengthfile][-lowmem][-m modelinitialfile][-many]"+
                                    "[-n numseq][-noautoopen][-nobed][-nobrowser][-nobrowserheader][-noenrich][-noimage][-nopseudolast][-p maxprocessors][-pseudo][-printposterior][-printstatebyline][-r maxiterations][-s seed][-scalebeta]"+
                                    "[-splitrows][-stateordering emission|transition]"+
                                   "[-t loadsmoothtransition][-u coorddir][-v anchorfiledir][-x maxseconds][-z zerotransitionpower] inputdir outputdir numstates assembly");
	    }
	} 
	else if (szcommand.equalsIgnoreCase("Reorder"))
        {
	    int ninitmethod = ChromHMM.INITMETHOD_INFORMATION;
	    String szoutfileID = "";
	    String szlabelmapping = null;
	    String szstateorderingfile = null;
	    String szcolumnorderingfile = null;
	    int nstateorder = ChromHMM.STATEORDER_FIXED;
	    boolean bnoordercols = false;
	    boolean bnoprintsegment = false;
	    boolean bprintimage = true;
	    int nargindex = 1;
	    String szreorderinbedfile = null;
            String szreorderoutbedfile = null;
	    boolean breordercolsmodel = false;

	    int nr=ChromHMM.DEFAULTCOLOR_R;
	    int ng=ChromHMM.DEFAULTCOLOR_G;
	    int nb=ChromHMM.DEFAULTCOLOR_B;

            try
	    {
	       while (nargindex < args.length-2)
	       {
                  if (args[nargindex].equals("-color"))
		  {
		     String szcolor = args[++nargindex];
		     StringTokenizer stcolor = new StringTokenizer(szcolor,",");
		     if (stcolor.countTokens()==3)
		     {
		        nr = Integer.parseInt(stcolor.nextToken().trim());
		        ng = Integer.parseInt(stcolor.nextToken().trim());
		        nb = Integer.parseInt(stcolor.nextToken().trim());
		     }
		     else
		     {
		        bok = false;
		     }
		  }
		  else if (args[nargindex].equals("-f"))
		  {
		      szcolumnorderingfile = args[++nargindex];
		  }
		  else if (args[nargindex].equals("-i"))
		  {
		     szoutfileID = args[++nargindex];
		  }
                  else if (args[nargindex].equals("-m"))
		  {
		      szlabelmapping = args[++nargindex];
		  }
                  else if (args[nargindex].equals("-noimage"))
		  {
		      bprintimage = false;
		  }
                  else if (args[nargindex].equals("-reordercolsmodelfile"))
		  {
		      breordercolsmodel = true;
		  }
		  else if (args[nargindex].equals("-o"))
		  {
		      szstateorderingfile = args[++nargindex];
		      nstateorder = ChromHMM.STATEORDER_USER;
		  }
                  else if (args[nargindex].equals("-r"))
	          {
		     szreorderinbedfile = args[++nargindex];
		     if (nargindex + 1 == args.length)
		     {
			 bok = false;
			 break;
		     }
		     else
		     {
		        szreorderoutbedfile = args[++nargindex];
		     }

		     if (szreorderinbedfile.equals(szreorderoutbedfile))
		     {
			 System.out.println("bedfilein cannot be the same as bedfileout");
			 bok = false;
		     }
		  }
		  else if (args[nargindex].equals("-stateordering"))
		  {
                     String sztoken = args[++nargindex];
                     if (sztoken.equals("emission"))
		     {
			 nstateorder = ChromHMM.STATEORDER_EMISSION;
		     } 
                     else if (sztoken.equals("transition"))
		     {
			 nstateorder = ChromHMM.STATEORDER_TRANSITION;
		     }
		     else
		     {
			 bok = false;
			 break;
		     }
		  }
		  else if (args[nargindex].equals("-holdcolumnorder"))
		  {
		     bnoordercols = true;
		  }
		  else
		  {
		     bok = false;
		  }
		  nargindex++;
	       }
	    }
	    catch (NumberFormatException ex)
	    {
		bok = false;
	    }

	    if ((szstateorderingfile == null) && (szreorderinbedfile != null))
	    {
	       bok = false;
	    }

	    if (bok&&(nargindex==args.length-2))
	    {
	       Color theColor = new Color(nr, ng, nb);
	       String szInitFile = args[nargindex++];
	       String szoutputdir = args[nargindex++];
	       File f = new File(szoutputdir);
	       if (!f.exists())
	       {
	          if (!f.mkdirs())
	          {
	             throw new IllegalArgumentException(szoutputdir+" does not exist and could not be created!");
	          }
	       }

	       if (bok)
	       {
	          boolean bprintsegments = !bnoprintsegment; 
	          boolean bordercols = !bnoordercols;
		  ChromHMM theHMM = new ChromHMM(szInitFile, szoutputdir, szstateorderingfile, szcolumnorderingfile, szoutfileID,nstateorder,
                                                 bordercols,theColor,szlabelmapping,bprintimage, szreorderinbedfile, 
                                                 szreorderoutbedfile, breordercolsmodel);
		  theHMM.reorderModel();
	       }
	    }
	    else
	    {
		bok = false;
	    }

	    if (!bok)
	    {

		System.out.println("usage: Reorder [-color r,g,b][-f columnorderingfile][-holdcolumnorder][-i outfileID]"+
                                   "[-m labelmappingfile][-noimage][-o stateorderingfile [-r bedfilein bedfileout]][-reordercolsmodelfile][-stateordering emission|transition] inputmodel outputdir");

	    }
	}
        else if (szcommand.equalsIgnoreCase("ConvertGeneTable"))
	{
            String path = ChromHMM.class.getProtectionDomain().getCodeSource().getLocation().getPath();
            String decodedPath = URLDecoder.decode(path, "UTF-8");
            String szprefixpath = decodedPath.substring(0, decodedPath.lastIndexOf("/") + 1);
            //String szfullprefixpathchromsizedir = szprefixpath+"/"+CHROMSIZESDIR;
            String szfullprefixpathcoorddir = szprefixpath+"/"+COORDDIR;
            String szfullprefixpathanchorfiledir = szprefixpath+"/"+ANCHORFILEDIR;
	    String szchromlengthfile = null;
            int nargindex = 1;
	    int npromoterwindow = 2000;
	    boolean bgzip = false;
	    boolean bnobin = false;
	    boolean bnoheader = false;
	    boolean bbiggenepred = false;
	    String szcoorddir;
	    String szanchordir;

            try
	    {
	       while (nargindex < args.length-3)
	       {
                  //else if (args[nargindex].equals("-u"))
		  //    {
		  //	  szfullprefixpathchromsizedir = args[++nargindex];
		  //    }
		  if (args[nargindex].equals("-biggenepred"))
		  {
		     bbiggenepred = true;
		  }
                  else if (args[nargindex].equals("-gzip"))
		  {
	       	     bgzip = true;
		  }
	          else if (args[nargindex].equals("-l"))
		  {
		     szchromlengthfile = args[++nargindex];
	          }
	          else if (args[nargindex].equals("-nobin"))
		  {
		      bnobin = true;
		  }
		  else if (args[nargindex].equals("-noheader"))
		  {
		      bnoheader = true;
		  }
                  else if (args[nargindex].equals("-u"))
		  {
		     szfullprefixpathcoorddir = args[++nargindex];
		  }
                  else if (args[nargindex].equals("-v"))
		  {
		     szfullprefixpathanchorfiledir = args[++nargindex];
	          }
	          else if (args[nargindex].equals("-w"))
		  {
		      npromoterwindow = Integer.parseInt(args[++nargindex]);
		  }
		  else
		  {
		     bok = false;
		  }
                  nargindex++;
	       }
	    }
            catch (NumberFormatException ex)
	    {
	       bok = false;
	    }

	    if ((bbiggenepred)&&(bnobin))
	    {
		bok = false;
	    }

            if (bok && (nargindex==args.length-3))
	    {
	       String sztable = args[nargindex++];
	       //String szchromlengths = args[nargindex++];
               String szprefix = args[nargindex++];
	       String szassembly = args[nargindex++];

	       szcoorddir = szfullprefixpathcoorddir + "/" + szassembly;
	       szanchordir = szfullprefixpathanchorfiledir + "/" + szassembly;

	       File f = new File(szcoorddir);
	       if (!f.exists())
	       {
	          if (!f.mkdirs())
	          {
	             throw new IllegalArgumentException(szcoorddir+" does not exist and could not be created!");
	          }
	       }

	       f = new File(szanchordir);
	       if (!f.exists())
	       {
	          if (!f.mkdirs())
	          {
	             throw new IllegalArgumentException(szanchordir+" does not exist and could not be created!");
	          }
	       }

	       if (szchromlengthfile == null)
	       {
		  //updated v.1.16
		  String szdefaultlengthfile = szprefixpath+"/"+CHROMSIZESDIR+"/"+szassembly+".txt";
		  File flength = new File(szdefaultlengthfile);
		  //File flength = new File(szfullprefixpathchromsizedir+"/"+szassembly+".txt");
		  if (flength.exists())
		  {
		      //szchromlengthfile = szfullprefixpathchromsizedir+"/"+szassembly+".txt";
		      szchromlengthfile = szdefaultlengthfile;
		  }
		  else
		  {
		      throw new IllegalArgumentException("chromosome length file not specified and default "+szdefaultlengthfile+" not found");
		  }
	       }


	       ConvertGeneTable.convertGeneTableToAnnotations(sztable, szprefix, szassembly, szcoorddir, szanchordir, szchromlengthfile, 
                                                              npromoterwindow, bgzip, bnobin, bnoheader,bbiggenepred);
	    }
	    else
	    {
		bok = false;
	    }

	    if (!bok)
	    {
                System.out.println("usage: ConvertGeneTable [-biggenepred|-nobin][-gzip][-l chromosomelengthfile][-noheader][-u coorddir][-v anchordir][-w promoterwindow] inputgenetable prefix assembly");
	    }

	}
	else
	{
	    System.out.println("Need to specify the mode BinarizeBam|BinarizeBed|BinarizeSignal|CompareModels|ConvertGeneTable|EvalSubset|LearnModel|MakeBrowserFiles"+
                               "|MakeSegmentation|MergeBinary|NeighborhoodEnrichment|StatePruning|OverlapEnrichment|Reorder|Version");

	}
    }
}

