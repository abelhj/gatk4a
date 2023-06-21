package  org.broadinstitute.hellbender.tools.walkers.haplotect;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.LocusWalker;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import htsjdk.samtools.util.Locatable;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.Comparator;
import java.util.LinkedHashSet;
import java.util.HashSet;
import java.util.Map;
import java.util.LinkedHashMap;
import java.util.ArrayList;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.File;
import java.lang.Runtime;

import org.broadinstitute.hellbender.tools.walkers.haplotect.haplotect_utils.*;


@CommandLineProgramProperties(
			      summary = "Generate contamination estimates from haplotype counts",
			      oneLineSummary = "Generate contamination estimates from haplotype counts",
			      programGroup = CoverageAnalysisProgramGroup.class)
@DocumentedFeature
public class Haplotect extends LocusWalker  {

    String version="0.3";

   
    @Argument(fullName = "minMappingQuality", shortName = "mmq", doc = "Minimum mapping quality of reads to count towards depth", optional = true, minValue = 0, maxValue = Integer.MAX_VALUE)
	int minMappingQuality = -1;
    @Argument(fullName = "minBaseQuality", shortName = "mbq", doc = "Minimum quality of bases to count towards depth", optional = true, minValue = 0, maxValue = Byte.MAX_VALUE)
	byte minBaseQuality = -1;
    @Argument(fullName = "pairsFile", shortName = "htp", doc = "file of locus pairs", optional = false)
	String pairsFile;
    @Argument(fullName = "debug", shortName = "debug", doc= "true for verbose output", optional = true)
	boolean debug = false;
    @Argument(fullName = "gstol", shortName= "gstol", doc = "error tolerance for golden section search", optional = true, minValue=0.00000001, maxValue=0.1)
	double gstol=0.005;
    @Argument(fullName = "outPrefix", shortName= "outPrefix", doc = "prefix for output files", optional = false)
	String outPrefix=null;
    @Argument(fullName = "uniform", shortName = "unif", doc= "use uniform prior for popn hap frequencies", optional = true)
	boolean uniform=false;
    @Argument(fullName = "minreads", shortName = "mr", doc = "minium number of reads at a locus to include", optional = true)
        int minreads=20;
    @Argument(fullName="minMultiPct", shortName = "minPct", doc="min pct multihaplotype sites", optional = true)
	int minMultiPct=5;
    @Argument(fullName="calcMLE", shortName = "calcmle", doc="false to skip mle calculation", optional = true)
        boolean calcMLE = true;
    
    ArrayList<SnpPair> pairs=null;
    LinkedHashSet<GenomeLoc> snps=null;    
    LinkedHashMap<SnpPair, LinkedHashMap<String, LinkedHashMap<Integer, BaseandQual> > > hapmap=null;
    LinkedHashMap<SnpPair, HapCounter> paircts=null;
    HashSet<SnpPair> multiPairs=null;

    GenomeLocParser gpl=null;
    PrintStream out=null;
    PrintStream log=null;
    String currentContig=null;
    int totalCounts=0;
    double contamCounts=0;
    int informativeSites=0;
    int totalSites=0;
    double aveContFrac=0;
    double meancov=0;
    

    @Override
    public void onTraversalStart() {
        BufferedReader br=null;
        pairs=new ArrayList<SnpPair>();
        snps=new LinkedHashSet<GenomeLoc>();
        gpl=new GenomeLocParser(this.getMasterSequenceDictionary());
        try {
            br=new BufferedReader(new FileReader(pairsFile));
            String curLine=null;
            while((curLine=br.readLine())!=null) {
                
                String[] spl=curLine.split("\\t");
                GenomeLoc gl1=gpl.parseGenomeLoc(spl[0]+":"+spl[1]);
                GenomeLoc gl2=gpl.parseGenomeLoc(spl[0]+":"+spl[2]);
		if( !(snps.contains(gl1) && snps.contains(gl2))) {
			snps.add(gl1);
			snps.add(gl2);
			pairs.add(new SnpPair(gl1, gl2, curLine, uniform));  //fix this !!
		} else {
		    System.err.println("Warning:  Duplicate snp pair "+gl1+" "+gl2+".  Using only first occurrence.");
		}
            }
	    br.close();
	    out=new PrintStream (new File(outPrefix+".haplotect.txt"));
	    log=new PrintStream (new File(outPrefix+".haplotectloci.txt"));
	    log.println("#Haplotect_v"+version); 
	    log.println("#"+pairs.size()+" SNP pairs read");
	    log.println("#chr\tSNP1\tSNP2\tall11\tall12\tall21\tall22\tpopn_counts\tdistance\ttotal_count\tsample_counts\t");
        }
        catch(IOException e) {
            e.printStackTrace();
        }
	
	currentContig="";

        hapmap=new LinkedHashMap<SnpPair, LinkedHashMap<String, LinkedHashMap<Integer, BaseandQual> > >();
	paircts=new LinkedHashMap<SnpPair, HapCounter>();
        multiPairs=new HashSet<SnpPair>();
        for(SnpPair pr : pairs) {
            hapmap.put(pr, new LinkedHashMap<String, LinkedHashMap<Integer, BaseandQual> >());
	}
    }

    @Override
    public void apply(AlignmentContext context, ReferenceContext ref, FeatureContext featureContext) {
        if (hasReference() && BaseUtils.isRegularBase(ref.getBase())) {
            Locatable gl1=context.getLocation();
	    GenomeLoc loc=gpl.parseGenomeLoc(gl1.getContig()+":"+gl1.getStart());
	    if(currentContig==null) {
	        currentContig=context.getContig();
            } else if(!context.getContig().equals(currentContig)) {
	        //System.err.println(loc+"\t"+Runtime.getRuntime().totalMemory()/(1024*1024)+"\t"+Runtime.getRuntime().freeMemory()/(1024*1024)+"\t"+Runtime.getRuntime().maxMemory()/(1024*1024));
                finishChrom();
	        currentContig=context.getContig();
	    }
            if(snps.contains(loc)) {
	        ReadPileup pileup = context.getBasePileup();
                LinkedHashMap<String, BaseandQual> readallele=new LinkedHashMap<String, BaseandQual>();
                for(PileupElement p : pileup) {
                    if(p.getMappingQual()> minMappingQuality && BaseUtils.isRegularBase(p.getBase()) && p.getQual()>minBaseQuality && !p.getRead().isDuplicate()) 
                        readallele.put(p.getRead().getName(), new BaseandQual(p));
                }
                LocusNameAllele value=new LocusNameAllele(loc, readallele);

                for(SnpPair pr: pairs) {
                    if(pr.matches(loc)) {
                        int order=pr.whichSnp(loc);
                        LinkedHashMap<String, BaseandQual> nts=value.getAlleles();
                        for(String readname : nts.keySet()) {
                            if(!hapmap.get(pr).containsKey(readname)) 
                                hapmap.get(pr).put(readname, new LinkedHashMap<Integer, BaseandQual>());
                            hapmap.get(pr).get(readname).put(order, nts.get(readname));
                        }    
                    }
                }
            }
	}
        return;
    }

  public void finishChrom() {
      System.err.println("finish chrom");
      if(currentContig==null) {
	  return;
      }
    for(SnpPair pr: hapmap.keySet()) {
      HapCounter hapcts=new HapCounter(pr);
      for(String name:hapmap.get(pr).keySet()) {
        if(hapmap.get(pr).get(name).keySet().size()==2) {               //only use reads covering both SNVs in pair
          hapcts.add(hapmap.get(pr).get(name).get(1), hapmap.get(pr).get(name).get(2));
          /*System.err.println("read\t"+name+"\t"+hapmap.get(pr).get(name).get(1)+"\t"+hapmap.get(pr).get(name).get(2));*/
        }
      }
      // requires minreads reads covering the locus to consider
      if (hapcts.totalCount()>minreads){		
        totalSites++;		
	paircts.put(pr, hapcts);
	meancov+=hapcts.totalCount();
	if(hapcts.numUniqueObsHaps()>2) {
	  multiPairs.add(pr);
          if(hapcts.fourthHapCount()>1) {
	      informativeSites++;
              contamCounts+=0.5*(hapcts.thirdHapCount()+hapcts.fourthHapCount());
              totalCounts+=hapcts.totalCount();
          } else if (hapcts.thirdHapCount()>1) {
              informativeSites++;
              contamCounts+=hapcts.thirdHapCount();
              totalCounts+=hapcts.totalCount();
          }
	  log.println(pr.pairInfo()+"\t"+pr.distance()+"\t"+hapcts.totalCount()+"\t"+hapcts.countString());
	} 
      }
      if(debug) {
	  hapcts.print(log);
	  
	  }
      hapmap.put(pr, new LinkedHashMap<String, LinkedHashMap<Integer, BaseandQual> >());
    }
  }

  @Override
  public Object onTraversalSuccess() {
        finishChrom();
	meancov/=totalSites;
	aveContFrac = 2.0*contamCounts/totalCounts;
	String[] spl=outPrefix.split("/");
	String id=spl[spl.length-1];
        out.println("#sample\ttotalSites\tinformativeSites\tmeanCoverage\tcontaminatingHaplotypeCount\ttotalReads\tcontaminationEstimate\tsnpEstimate\tsnpEstimateCI");
        if (!calcMLE) {
	    out.println(id+"\t"+totalSites+"\t"+informativeSites+"\t"+String.format("%.3f",meancov)+"\t"+String.format("%.1f", contamCounts)+"\t"+totalCounts+"\t"+String.format("%.4f", aveContFrac)+"\tNA\tNA");
        } else {
	  System.err.println("Calcuating mle");
	  PairLogLik mleCalculator=new PairLogLik(paircts, gstol, debug, log);
	  System.err.println("Calcuating mle");
	  mleCalculator.calcMleCI();	
	  double mle=mleCalculator.getMle();
	  double[] CI=mleCalculator.getCI();
	  out.println(id+"\t"+totalSites+"\t"+informativeSites+"\t"+String.format("%.3f",meancov)+"\t"+String.format("%.1f", contamCounts)+"\t"+totalCounts+"\t"+String.format("%.4f", aveContFrac)+"\t"+String.format("%.3f",mle)+"\t"+String.format("%.3f", CI[0])+"-"+String.format("%.3f", CI[1]));
        }
        return "";
  }
}
