package org.broadinstitute.hellbender.tools.walkers.haplotect.haplotect_utils;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;

public class BaseandQual{

    private byte base;
    private byte qual;
   
    public BaseandQual(PileupElement pel) {
	base=pel.getBase();
	qual=pel.getQual();
    }

    public byte getBase() {
	return base;
    }

    public byte getQual() {
	return qual;
    }

    public String toString() {
	return ""+base+"_"+qual;
    }
}