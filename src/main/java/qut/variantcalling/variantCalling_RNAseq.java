/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Project/Maven2/JavaApp/src/main/java/${packagePath}/${mainClassName}.java to edit this template
 */

package qut.variantcalling;
import htsjdk.samtools.*;
import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import picard.sam.CreateSequenceDictionary;
import qut.variantcalling.including.Nucl_str;

/**
 *
 * @author an
 */
public class variantCalling_RNAseq {

    public static void main(String[] args) throws IOException {
        int min_total_read_depth = 100;
        String gff3_fn = "";
        String chromosome_targeted = "all";
        String bam_fn = "", genome_fn = "";

        if (System.getProperty("os.name").startsWith("Windows")) {
            bam_fn = "C:\\Jiyuan\\sourceCode\\Zuba\\NbLab360.genome.fastaOZBenth_all_R1.fastp.fq.gz.sorted.bam_NbLab360C01.bam";
            bam_fn = "C:\\Jiyuan\\sourceCode\\RNASeqBrowser\\data\\NbLab360\\t.bam";
            bam_fn = "C:\\Jiyuan\\sourceCode\\RNASeqBrowser\\data\\NbLab360\\QLD_read_LAB360_genome.chr01_partial.bam";
            bam_fn = "C:\\Jiyuan\\sourceCode\\RNASeqBrowser\\data\\NbLab360\\Lall.bwa.Lab360genome.nochrom0.bam.sorted.bam.NbLab360C01.bam";
            bam_fn = "C:\\Jiyuan\\sourceCode\\zuba\\LAB.genome.Lall.Aligned.sortedByCoord.out.bam.NbLab360C01.bam";
            bam_fn = "C:\\Jiyuan\\sourceCode\\RNASeqBrowser\\data\\NbQld183\\Qall.bwa.QLD183genome.nochrom0.bam.sorted.bam.NbQld183C01.bam";
            bam_fn = "C:\\Jiyuan\\sourceCode\\zuba\\NbQld183C15.10M.bam";
            bam_fn = "C:\\Jiyuan\\sourceCode\\zuba\\NbQld183C10.12m.bam";
            bam_fn = "C:\\Jiyuan\\sourceCode\\zuba\\NbQld183C10.10M.bam";
            bam_fn = "C:\\Jiyuan\\sourceCode\\maven\\vaiantCalling\\demoData\\NbQld183C10.1M.bam";
            genome_fn = "C:\\Jiyuan\\sourceCode\\NB_annotation\\LAB360\\NbLab360.genome.fasta";
            genome_fn = "C:\\Jiyuan\\sourceCode\\NB_annotation\\QLD183\\NbQld183.genome.fasta";
            genome_fn = "C:\\Jiyuan\\sourceCode\\zuba\\C10.NbQld183.fasta";
            genome_fn = "C:\\Jiyuan\\sourceCode\\maven\\vaiantCalling\\demoData\\C10.NbQld183.fasta";
            gff3_fn = "C:\\Jiyuan\\sourceCode\\NB_annotation\\LAB360\\NbLab360.hairTail.gff3.addAnnot_blast.gff3";
            gff3_fn = "C:\\Jiyuan\\sourceCode\\NB_annotation\\QLD183\\NbQld183.v103.gff3";
            gff3_fn = "C:\\Jiyuan\\sourceCode\\zuba\\NbQld183.v103.gff3";
            gff3_fn = "C:\\Jiyuan\\sourceCode\\maven\\vaiantCalling\\demoData\\NbQld183.v103.gff3";
//            chromosome_targeted = "NbQld183C15";
        } else {
            if (args.length < 2) {
                error_commandline();
            }
            for (int i = 0; i < args.length - 2; i++) {
                if (args[i].equals("--gff3") || args[i].equals("-g")) {
                    gff3_fn = args[++i];
                }else if (args[i].equals("--chr-target") || args[i].equals("-c")) {
                    chromosome_targeted = args[++i];
                }else if (args[i].equals("--min_total_read_depth") || args[i].equals("-d")) {
                    min_total_read_depth = Integer.parseInt(args[++i]);
                }else{
                    error_commandline();
                }
            }
            bam_fn = args[args.length - 2];
            genome_fn = args[args.length - 1];
        }
        new variantCalling_RNAseq().proc(bam_fn, genome_fn, chromosome_targeted, min_total_read_depth, gff3_fn);
    }

    private static void error_commandline() {
        System.err.println("java -cp vaiantCalling.jar qut.vaiantcalling.variantCalling_RNAseq \\");
        System.err.println("--gff3/-g xxx.gff3 \\");
        System.err.println("--chr-target/-c chrX \\");
        System.err.println("--min_total_read_depth/-d 100 \\");
        System.err.println("bam_filename genome_filename");
        System.exit(0);
    }
    
    void createDict(String genome_fn){
        
        // Create the list of arguments for the Picard tool
        List<String> picardArgs = new ArrayList<>();
        picardArgs.add("R=" + genome_fn);  // Input reference genome file
        picardArgs.add("O=" + genome_fn.replace("fasta", "dict"));       // Output dictionary file

        // Create an instance of CreateSequenceDictionary and run it with the arguments
        CreateSequenceDictionary createSequenceDictionary = new CreateSequenceDictionary();
        createSequenceDictionary.instanceMain(picardArgs.toArray(new String[0]));

    }
    
    void proc(String bam_fn, String genome_fn, String chromosome_targeted, int min_total_read_depth, String gff3_fn) throws IOException{
        String pure_fn = genome_fn.substring(genome_fn.replace("\\", "/").lastIndexOf("/")+1);
        String middle_fn = bam_fn + "."+pure_fn+".varant."+chromosome_targeted+".txt";
        BufferedWriter bw = new BufferedWriter(new FileWriter(middle_fn));
        bw.write("chromosome\t" + "position(o-based)\t" + "ref_bp:#reads\t" + "SNP_bp:#reads;...\t"+ "Insert_bp:#reads;...\t" + "Del_bp:#reads;...\n");
        BufferedWriter bw_detail = new BufferedWriter(new FileWriter(middle_fn+".detail.txt"));
        bw_detail.write("readId\t" + "chromosome\t" + "position(o-based)\t" + "ref_bp\t" + "read_bp\n");
        SamReader reader = SamReaderFactory.makeDefault().open(new File(bam_fn));
        if (!new File(genome_fn.replace(".fasta", ".dict")).exists()) {
            System.err.println("please create genome sequence index using picard and samtools");
        }
        ReferenceSequenceFile refFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(new File(genome_fn));
        Map<String, Integer> chrom_size = new TreeMap<>();
        for (SAMSequenceRecord sequenceRecord : refFile.getSequenceDictionary().getSequences()) {
            String sequenceName = sequenceRecord.getSequenceName();
            int sequenceLength = sequenceRecord.getSequenceLength();
            System.out.println(sequenceName + "\t" + sequenceLength);
            chrom_size.put(sequenceName, sequenceLength);
        }
        String current_chr = "";
        Nucl_str[] Nucl_current_chr = new Nucl_str[0];
        for (SAMRecord read : reader) {
            if(!chromosome_targeted.equals("all") && !read.getReferenceName().contains(chromosome_targeted)){
                continue;
            }
            if (!read.getReferenceName().equals(current_chr)) {
                for (int i = 0; i < Nucl_current_chr.length; i++) {
                    Nucl_str Nucl_current_chr1 = Nucl_current_chr[i];
                    if (Nucl_current_chr1.alter != null) {
                        bw.write(current_chr + "\t" + i + "\t" + Nucl_current_chr1.ref + ":" + Nucl_current_chr1.count);
                        String delimer = "\t";
                        for (char alt : Nucl_current_chr1.alter.keySet()) {
                            bw.write(delimer + alt + ":" + Nucl_current_chr1.alter.get(alt));
                            delimer = ";";
                        }
                        if (Nucl_current_chr1.Ins != null) {
                            delimer = "\t";
                            for (String ins : Nucl_current_chr1.Ins.keySet()) {
                                bw.write(delimer + ins + ":" + Nucl_current_chr1.Ins.get(ins));
                                delimer = ";";
                            }
                        } else {
                            bw.write("\tN/A");
                        }
                        if (Nucl_current_chr1.Del != null) {
                            delimer = "\t";
                            for (String del : Nucl_current_chr1.Del.keySet()) {
                                bw.write(delimer + del + ":" + Nucl_current_chr1.Del.get(del));
                                delimer = ";";
                            }
                        } else {
                            bw.write("\tN/A");
                        }
                        bw.write("\n");
                    }
                }
                
                System.out.println(current_chr);
                if(chrom_size.get(read.getReferenceName()) == null){
                    System.err.println("chromosome:"+read.getReferenceName() + " is not exist");
                    continue;
                }
                current_chr = read.getReferenceName();
                Nucl_current_chr = new Nucl_str[chrom_size.get(current_chr)];
                String refSeq = refFile.getSubsequenceAt(current_chr, 1, chrom_size.get(current_chr)).getBaseString();
                for(int i = 0; i < refSeq.length(); i++){
                    Nucl_current_chr[i] = new Nucl_str(refSeq.charAt(i));
                }
            }
            if (!read.getReadUnmappedFlag() && !read.isSecondaryOrSupplementary()) {
                Cigar cigar = read.getCigar();
                int readPos = 0;
                int refPos = read.getAlignmentStart() - 1;
                for (CigarElement cigarElement : cigar.getCigarElements()) {
                    // Check if the operation is a SNP
                    switch (cigarElement.getOperator()) {
                        case M:
                            ReferenceSequence refSeq = refFile.getSubsequenceAt(
                                    read.getReferenceName(), refPos + 1, refPos + cigarElement.getLength());
                            byte[] refBases = refSeq.getBases();
                            byte[] readBases = read.getReadBases();
                            for (int i = 0; i < cigarElement.getLength() - 0; i++) {//the two side of read not considering
                                    if (refBases[i] != readBases[readPos + i]) {
                                        bw_detail.write(read.getReadName() + "\t" + read.getReferenceName() + "\t"
                                                + (refPos + i) + "\t"
                                                + (char) (refBases[i]) + "\t"
                                                + (char) (readBases[readPos + i]) + "\n");
                                        char c = (char) readBases[readPos + i];
                                        if (Nucl_current_chr[refPos + i].alter == null) {
                                            Nucl_current_chr[refPos + i].alter = new TreeMap<>();
                                        }
                                        if (Nucl_current_chr[refPos + i].alter.get(c) == null) {
                                            Nucl_current_chr[refPos + i].alter.put(c, 1);
                                        } else {
                                            Nucl_current_chr[refPos + i].alter.put(c, Nucl_current_chr[refPos + i].alter.get(c) + 1);
                                        }
                                    } else {
                                        Nucl_current_chr[refPos + i].count++;
                                    }
                                
                            }
                            readPos += cigarElement.getLength();
                            refPos += cigarElement.getLength();
                            break;
                        case N:
                            refPos += cigarElement.getLength();
                            break;
                        case D:
                            ReferenceSequence refSeq_del = refFile.getSubsequenceAt(
                                    read.getReferenceName(), refPos + 1, refPos + cigarElement.getLength());
                            String del = refSeq_del.getBaseString();
                            if (Nucl_current_chr[refPos].Del == null) {
                                Nucl_current_chr[refPos].Del = new TreeMap<>();
                            }
                            if (Nucl_current_chr[refPos].Del.get(del) == null) {
                                Nucl_current_chr[refPos].Del.put(del, 1);
                            } else {
                                Nucl_current_chr[refPos].Del.put(del, Nucl_current_chr[refPos].Del.get(del) + 1);
                            }
                            refPos += cigarElement.getLength();
                            break;
                        case I:
                            String ins = new String(Arrays.copyOfRange(read.getReadBases(), readPos, readPos + cigarElement.getLength()));
                            if (Nucl_current_chr[refPos].Ins == null) {
                                Nucl_current_chr[refPos].Ins = new TreeMap<>();
                            }
                            if (Nucl_current_chr[refPos].Ins.get(ins) == null) {
                                Nucl_current_chr[refPos].Ins.put(ins, 1);
                            } else {
                                Nucl_current_chr[refPos].Ins.put(ins, Nucl_current_chr[refPos].Ins.get(ins) + 1);
                            }
                            readPos += cigarElement.getLength();
                            break;
                        case S:
                            readPos += cigarElement.getLength();
                            break;
                        case H:
                            // do nothing since hard clip does not affect the read or reference position
                            break;
                        default:
                            // handle other operators as needed
                            break;
                    }
                }
            }
        }
        for (int i = 0; i < Nucl_current_chr.length; i++) {
            Nucl_str Nucl_current_chr1 = Nucl_current_chr[i];
            if (Nucl_current_chr1.alter != null) {
                bw.write(current_chr+"\t"+i + "\t" + Nucl_current_chr1.ref + ":" + Nucl_current_chr1.count);
                String delimer = "\t";
                for (char alt : Nucl_current_chr1.alter.keySet()) {
                    bw.write(delimer + alt + ":" + Nucl_current_chr1.alter.get(alt));
                    delimer = ";";
                }
                if (Nucl_current_chr1.Ins != null) {
                    delimer = "\t";
                    for (String ins : Nucl_current_chr1.Ins.keySet()) {
                        bw.write(delimer + ins + ":" + Nucl_current_chr1.Ins.get(ins));
                        delimer = ";";
                    }
                }else{
                    bw.write("\tN/A");
                }
                if (Nucl_current_chr1.Del != null) {
                    delimer = "\t";
                    for (String del : Nucl_current_chr1.Del.keySet()) {
                        bw.write(delimer + del + ":" + Nucl_current_chr1.Del.get(del));
                        delimer = ";";
                    }
                }else{
                    bw.write("\tN/A");
                }
                bw.write("\n");
            }
        }
        reader.close();
        refFile.close();
        bw.close();
        bw_detail.close();
        new postVariantCalling(middle_fn, gff3_fn, min_total_read_depth).proc();
    }
}
