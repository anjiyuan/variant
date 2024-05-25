/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package qut.variantcalling.including;

/**
 *
 * @author an
 */
public class loci {
    public String chr;
    public char strand;
    public int start;
    public int end;
    public StringBuilder seq;

    public loci(String chr, char strand, int start, int end) {
        this.chr = chr;
        this.strand = strand;
        this.start = start;
        this.end = end;
    }
    
    public loci(loci l) {
        this.chr = l.chr;
        this.strand = l.strand;
        this.start = l.start;
        this.end = l.end;
    }

    public int len(){
        return end - start;
    }
    @Override
    public String toString() {
        return "loci{" + "chr=" + chr + ", strand=" + strand + ", start=" + start + ", end=" + end + '}';
    }
    
}
