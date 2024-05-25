/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package qut.variantcalling.including;


import java.util.Map;

/**
 *
 * @author an
 */
public class Nucl_str {
    public char ref;
    public int count;
    public Map<Character, Integer> alter;
    public Map<String, Integer> Ins;
    public Map<String, Integer> Del;

    public Nucl_str(char Nt) {
        this.ref = Nt;
        this.count = 0;
    }

    @Override
    public String toString() {
        return "Nucl_str{" + "ref=" + ref + ", count=" + count + ", alter=" + alter + '}';
    }

    public int alt_count(){
        int ret = 0;
        if (alter != null) {
            for (char key : alter.keySet()) {
                ret += alter.get(key);
            }
        }
        return ret;
    }
}
