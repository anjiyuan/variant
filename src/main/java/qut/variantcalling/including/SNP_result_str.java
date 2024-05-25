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
public class SNP_result_str {
    public String chr;
    public int pos;
    public String ref;
    public String alt;
    public String Ins;
    public String Del;

    public SNP_result_str(String line) {
        String[] strarray = line.split("\t");
        if (strarray.length == 6) {
            this.chr = strarray[0];
            this.pos = Integer.parseInt(strarray[1]);
            this.ref = strarray[2];
            this.alt = strarray[3];
            this.Ins = strarray[4];
            this.Del = strarray[5];
        }
    }
    public int total_depth(){
        int ret = Integer.parseInt(ref.split(":")[1]);
        for(String str : alt.split(";")){
            ret += Integer.parseInt(str.split(":")[1]);
        }
        if(!Ins.equals("N/A"))
        for(String str : Ins.split(";")){
            ret += Integer.parseInt(str.split(":")[1]);
        }
        if(!Del.equals("N/A"))
        for(String str : Del.split(";")){
            ret += Integer.parseInt(str.split(":")[1]);
        }
        return ret;
    }
    public int count_ref(){
        return Integer.parseInt(ref.split(":")[1]);
    }

    public String count_alt() {
        int max_count = 0;
        char Nuc_max_count = ' ';
        for (String str : alt.split(";")) {
            if (max_count < Integer.parseInt(str.split(":")[1])) {
                max_count = Integer.parseInt(str.split(":")[1]);
                Nuc_max_count = str.split(":")[0].charAt(0);
            }
        }
        return Nuc_max_count + ":" + max_count;
    }
}
