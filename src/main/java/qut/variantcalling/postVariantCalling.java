/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package qut.variantcalling;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.TreeMap;
import qut.variantcalling.including.SNP_result_str;
import qut.variantcalling.including.loci;

/**
 *
 * @author an
 */
public class postVariantCalling {
    int min_total_read_depth = 100;
    String fn = "";
    String gff3_fn;
    public static void main(String[] args) throws IOException {
        String fn1;
        String gff3_fn1 = "";
        int min_deep1 = 100;
        if (System.getProperty("os.name").startsWith("Windows")) {
            fn1 = "C:\\Jiyuan\\sourceCode\\zuba\\Qall.bwa.QLD183genome.nochrom0.bam.sorted.bam.NbQld183.genome.fasta.extractSNP.07.C01.txt";
            fn1 = "C:\\Jiyuan\\sourceCode\\RNASeqBrowser\\data\\NbLab360\\Lall.bwa.Lab360genome.nochrom0.bam.sorted.bam.NbLab360.genome.fasta.extractSNP.07.C01.txt";
            fn1 = "C:\\Jiyuan\\sourceCode\\RNASeqBrowser\\data\\NbLab360\\Lall.bwa.Lab360genome.nochrom0.bam.sorted.bam.NbLab360C01.bam.NbLab360.genome.fasta.extractSNP.08.all.txt";
            fn1 = "C:\\Jiyuan\\sourceCode\\zuba\\QLD.genome.Lall.Aligned.sortedByCoord.out.bam.NbQld183.genome.fasta.extractSNP.08.C01.txt";
            fn1 = "C:\\Jiyuan\\sourceCode\\zuba\\OZBenth2_all_R1.fastp.fq.gz.Lab360_genome.nochr00.bam.sorted.bam.NbLab360.genome.fasta.extractSNP.08.C01.txt";
            gff3_fn1 = "C:\\Jiyuan\\sourceCode\\RNASeqBrowser\\data\\NbLab360\\NbLab360.v103.gff3";
            gff3_fn1 = "C:\\Jiyuan\\sourceCode\\NB_annotation\\LAB360\\NbLab360.hairTail.gff3.addAnnot_blast.gff3";
            gff3_fn1 = "C:\\Jiyuan\\sourceCode\\NB_annotation\\QLD183\\NbQld183.hairTail.gff3.addAnnot_blast.gff3";
            gff3_fn1 = "";
        } else {
            if (args.length < 1) {
                error_commandline();
            }
            for (int i = 0; i < args.length - 2; i++) {
                if (args[i].equals("--gff3") || args[i].equals("-g")) {
                    gff3_fn1 = args[++i];
                }
                if (args[i].equals("--min_total_read_depth") || args[i].equals("-d")) {
                    min_deep1 = Integer.parseInt(args[++i]);
                }
            }
            fn1 = args[args.length - 1];
        }
        new postVariantCalling(fn1, gff3_fn1, min_deep1).proc();
    }

    private static void error_commandline() {
        System.err.println("java -cp vaiantCalling.jar qut.vaiantcalling.variantCalling_RNAseq \\");
        System.err.println("--gff3/-g xxx.gff3 \\");
        System.err.println("--min_total_read_depth/-d 100 \\");
        System.err.println("variantCalling result file");
        System.exit(0);
    }

    public postVariantCalling(String fn,String gff3_fn, int min_total_read_depth) {
        this.fn = fn;
        this.gff3_fn = gff3_fn;
        this.min_total_read_depth = min_total_read_depth;
    }

    public postVariantCalling(String fn) {
        this.fn = fn;
        this.min_total_read_depth = 100;
    }

    boolean InExon_intron(List<loci> ls, int pos) {
        for (loci ll : ls) {
            if ((ll.end > pos) && (ll.start < pos)) {
                return true;
            }
        }
        return false;
    }

    public void proc() throws IOException {
        TreeMap<String, List<loci>> exons = new TreeMap<>();
        TreeMap<String, List<loci>> introns = new TreeMap<>();
        if (!gff3_fn.isEmpty()) {
            BufferedReader gff3 = new BufferedReader(new FileReader(gff3_fn));
            String l;
            int pre_intron_end_or_start = -1;
            while ((l = gff3.readLine()) != null) {
                String[] strarray = l.split("\t");
                if ((strarray.length > 2) && strarray[2].equals("exon")) {
                    List<loci> tmp1 = exons.get(strarray[0]);
                    if (tmp1 == null) {
                        tmp1 = new ArrayList<>();
                    }
                    tmp1.add(new loci(strarray[0], strarray[6].charAt(0), Integer.parseInt(strarray[3]) - 1, Integer.parseInt(strarray[4])));
                    exons.put(strarray[0], tmp1);
                    if (pre_intron_end_or_start > -1) {//not first exon
                        List<loci> tmp = introns.get(strarray[0]);
                        if (tmp == null) {
                            tmp = new ArrayList<>();
                        }
                        if (strarray[6].equals("+")) {
                            tmp.add(new loci(strarray[0], strarray[6].charAt(0), pre_intron_end_or_start, Integer.parseInt(strarray[3])));
                        } else {
                            tmp.add(new loci(strarray[0], strarray[6].charAt(0), Integer.parseInt(strarray[4]), pre_intron_end_or_start));
                        }
                        introns.put(strarray[0], tmp);
                    }
                    if (strarray[6].equals("+")) {
                        pre_intron_end_or_start = Integer.parseInt(strarray[4]);
                    } else {
                        pre_intron_end_or_start = Integer.parseInt(strarray[3]);
                    }
                } else if ((strarray.length > 2) && (strarray[2].equals("gene") || strarray[2].equals("mRNA") || strarray[2].equals("transcript"))) {
                    pre_intron_end_or_start = -1;
                }
            }

            gff3.close();
        }
        
        BufferedReader br = new BufferedReader(new FileReader(fn));
//        BufferedWriter bw_summary = new BufferedWriter(new FileWriter(fn+".summary"));
        BufferedWriter bw_mismatch = new BufferedWriter(new FileWriter(fn + ".homozygous"));
        bw_mismatch.write("chromosome\t" + "position(o-based)\t" + "ref_bp:#reads\t" + "SNP_bp:#reads;...\t" + "Insert_bp:#reads;...\t" + "Del_bp:#reads;...\n");
        BufferedWriter bw_heterozygous = new BufferedWriter(new FileWriter(fn + ".heterozygous"));
        bw_heterozygous.write("chromosome\t" + "position(o-based)\t" + "ref_bp:#reads\t" + "SNP_bp:#reads;...\t" + "Insert_bp:#reads;...\t" + "Del_bp:#reads;...\n");
        String line;
        String header = br.readLine();
        while ((line = br.readLine()) != null){
            SNP_result_str srs = new SNP_result_str(line);
//            bw_summary.write(line + "\n");
            if((srs.total_depth() >= min_total_read_depth)){
                if (!gff3_fn.isEmpty()) {
                    if (InExon_intron(exons.get(srs.chr), srs.pos) && !InExon_intron(introns.get(srs.chr), srs.pos)) {
                        if (srs.count_ref() < srs.total_depth() * 0.2) {
                            bw_mismatch.write(line + "\n");
                        } else if (srs.count_ref() < srs.total_depth() * 0.8) {
                            bw_heterozygous.write(line + "\n");
                        }
                    }
                } else {
                    if (srs.count_ref() < srs.total_depth() * 0.2) {
                        bw_mismatch.write(line + "\n");
                    } else if (srs.count_ref() < srs.total_depth() * 0.8) {
                        bw_heterozygous.write(line + "\n");
                    }
                }
            }
        }
        br.close();
//        bw_summary.close();
        bw_mismatch.close();
        bw_heterozygous.close();
    }
}
