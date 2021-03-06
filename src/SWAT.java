import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class SWAT {

	// Defaults [killed by main]

	private int gop = -11;
	private int gep = -1;
	private int sfsp = -20;
	private int dfsp = -40;
	private int stp = -4;
	private String matrix = "blosum62";
	private int[][] sc_matrix;
	private String name1 = "sequence 1(p)";
	private String name2 = "sequence 2(n)";
	private String special = "";
	private boolean wc_b_a = false;

	static int res = 78;

	// Defaults end

	private String sequence1;
	private String sequence2;
	private String seq1_al = "";
	private String alig_al = "";
	private String alig_o = "";
	private String seq2o_al = "";
	private String seq2_al = "";
	private int len1, len2, len_al;
	private int[][] s;
	private int[][] e;
	private int[][] h;
	private Cell[][] c_way;
	private Mutator mut;
	private int[] wb_counter;
	private int bestrow, bestcol;
	private int bestval;
	private int startrow, startcol;

	public SWAT() {
	}

	public SWAT(String n1, String seq1, String n2, String seq2, int gop, int gep, int sfsp, int dfsp, int stp, String matrix, String special, boolean worst) throws IOException {
		this.name1 = n1;
		this.sequence1 = "|||" + seq1.toUpperCase().replaceAll(" ", "");
		this.name2 = n2;
		this.sequence2 = "|||" + seq2.toUpperCase().replaceAll(" ", "");
		this.len1 = sequence1.length();
		this.len2 = sequence2.length();
		this.gop = gop;
		this.gep = gep;
		this.sfsp = sfsp;
		this.dfsp = dfsp;
		this.matrix = matrix;
		// this.sc_matrix = Matrix.mgen("./matrices/" + matrix);
		this.sc_matrix = Matrix.mgen(matrix, stp);
		if (special.equals("")) {
			this.special = "-";
		} else {
			this.special = special;
		}
		this.wc_b_a = worst;
	}

	// nuc-triplet 2 Aminosaeure
	private static String convertN2A(String nuc) {
		nuc = nuc.toUpperCase().replaceAll("T", "U");
		String p = "";
		for (int i = 0; i < nuc.length() - 3; i += 3)
			p = p + Allookup.trilookup(nuc.substring(i, i + 3));

		return p;

	}

	public void align_wb() {
		String p = sequence1;
		// String n = sequence2;
		String n = sequence2.replaceAll("T", "U");

		s = new int[len1 + 1][len2 + 1];
		e = new int[len1 + 1][len2 + 1];
		h = new int[len1 + 1][len2 + 1];
		c_way = new Cell[len1 + 1][len2 + 1];
		Cell c1 = new Cell();
		Cell c2 = new Cell();
		Cell c3 = new Cell();
		Cell c4 = new Cell();
		Cell c5 = new Cell();
		for (int i = 0; i <= len1; i++) {
			c_way[i][0] = new Cell();
			s[i][0] = 0;
			e[i][0] = 0;
			h[i][0] = 0;
		}
		for (int j = 0; j <= len2; j++) {
			c_way[0][j] = new Cell();
			s[0][j] = 0;
			e[0][j] = 0;
			h[0][j] = 0;
		}

		for (int i = 1; i < len1; i++) {
			for (int j = 5; j < len2; j++) {
				e[i][j] = Math.max(s[i - 1][j] + gop + gep, e[i - 1][j] + gep);
				h[i][j] = Math.max(s[i][j - 3] + gop + gep, h[i][j - 3] + gep);

				// Check auf wild bases:
				// REG-EXPRESSION BASTELN FÜR CHECK, kleine Optimierung für
				// Algo!

				// wb at c1 (Speichern, was wildbase war!!! in
				// original-Sequence!!!)
				if (wb_check("" + n.charAt(j))) {
					c1 = wbm1(i - 1, j - 1, p.charAt(i), n.charAt(j));
				} else {
					c1 = m1(i - 1, j - 1, p.charAt(i), n.charAt(j));
				}
				// wb at c2
				if (wb_check("" + n.charAt(j - 1) + n.charAt(j))) {
					c2 = wbm2(i - 1, j - 1, p.charAt(i), "" + n.charAt(j - 1) + n.charAt(j));
				} else {
					c2 = m2(i - 1, j - 2, p.charAt(i), "" + n.charAt(j - 1) + n.charAt(j));
				}
				// wb at c3
				if (wb_check("" + n.charAt(j - 2) + n.charAt(j - 1) + n.charAt(j))) {
					c3 = wbm3(i - 1, j - 3, p.charAt(i), "" + n.charAt(j - 2) + n.charAt(j - 1) + n.charAt(j));
					// c3.nucorg = "" + n.charAt(j - 2) + n.charAt(j - 1) +
					// n.charAt(j); // test für wb am Ende
				} else {
					c3 = m3(i - 1, j - 3, p.charAt(i), "" + n.charAt(j - 2) + n.charAt(j - 1) + n.charAt(j));
				}
				// wb at c4
				if (wb_check("" + n.charAt(j - 3) + n.charAt(j - 2) + n.charAt(j - 1) + n.charAt(j))) {
					c4 = wbm4(i - 1, j - 4, p.charAt(i), "" + n.charAt(j - 3) + n.charAt(j - 2) + n.charAt(j - 1) + n.charAt(j));
				} else {
					c4 = m4(i - 1, j - 4, p.charAt(i), "" + n.charAt(j - 3) + n.charAt(j - 2) + n.charAt(j - 1) + n.charAt(j));
				}
				// wb at c5
				if (wb_check("" + n.charAt(j - 4) + n.charAt(j - 3) + n.charAt(j - 2) + n.charAt(j - 1) + n.charAt(j))) {
					c5 = wbm5(i - 1, j - 5, p.charAt(i), "" + n.charAt(j - 4) + n.charAt(j - 3) + n.charAt(j - 2) + n.charAt(j - 1) + n.charAt(j));
				} else {
					c5 = m5(i - 1, j - 5, p.charAt(i), "" + n.charAt(j - 4) + n.charAt(j - 3) + n.charAt(j - 2) + n.charAt(j - 1) + n.charAt(j));
				}

				// c4 = m4(i - 1, j - 4, p.charAt(i), "" + n.charAt(j - 3) +
				// n.charAt(j - 2) + n.charAt(j - 1) + n.charAt(j));
				// c5 = m5(i - 1, j - 5, p.charAt(i), "" + n.charAt(j - 4) +
				// n.charAt(j - 3) + n.charAt(j - 2) + n.charAt(j - 1) +
				// n.charAt(j));

				// check best way start
				s[i][j] = Math.max(Util.max(s[i - 1][j - 1] + c1.score, s[i - 1][j - 2] + c2.score, s[i - 1][j - 3] + c3.score, s[i - 1][j - 4] + c4.score),
						Util.max(s[i - 1][j - 5] + c5.score, e[i][j], h[i][j], 0));

				if (s[i][j] == s[i - 1][j - 1] + c1.score) {
					c_way[i][j] = c1;
				}
				if (s[i][j] == s[i - 1][j - 2] + c2.score) {
					c_way[i][j] = c2;
				}
				if (s[i][j] == s[i - 1][j - 3] + c3.score) {
					c_way[i][j] = c3;
				}
				if (s[i][j] == s[i - 1][j - 4] + c4.score) {
					c_way[i][j] = c4;
				}
				if (s[i][j] == s[i - 1][j - 5] + c5.score) {
					c_way[i][j] = c5;
				}

				if (s[i][j] == 0) {
					c_way[i][j] = new Cell();
					c_way[i][j].casus = -42;
				}
				if (s[i][j] == e[i][j]) {
					c_way[i][j] = new Cell();
					c_way[i][j].casus = 7;
					c_way[i][j].p_row = i - 1;
					c_way[i][j].p_col = j;
					c_way[i][j].acid = p.charAt(i);
					c_way[i][j].nuc = "---";
					c_way[i][j].triplet = "---";
				}
				if (s[i][j] == h[i][j]) {
					c_way[i][j] = new Cell();
					c_way[i][j].casus = 8;
					c_way[i][j].p_row = i;
					c_way[i][j].p_col = j - 3;
					c_way[i][j].acid = '_';
					c_way[i][j].nuc = "" + n.charAt(j - 2) + n.charAt(j - 1) + n.charAt(j);
					c_way[i][j].triplet = "" + n.charAt(j - 2) + n.charAt(j - 1) + n.charAt(j);
				}
				if (s[i][j] > bestval) {
					bestval = s[i][j];
					bestrow = i; // row;
					bestcol = j;
				}
				// check best way end
			}
		}
		// System.out.println("alignment score:" + bestval);
	}

	public Cell wbm1(int pi, int pj, char a, char n) {
		Cell wbx = m1(pi, pj, a, Allookup.wb_lookup(n).charAt(0));
		Cell wby = new Cell();
		for (int q = 1; q < Allookup.wb_lookup(n).length(); q++) {
			wby = m1(pi, pj, a, Allookup.wb_lookup(n).charAt(q));
			if (mmMut(wby.score, wbx.score, this.wc_b_a)) {
				wbx = wby;
			}
		}
		switch (wbx.casus) {
		case 11:
			wbx.nucorg = n + "  ";
			break;
		case 12:
			wbx.nucorg = " " + n + " ";
			break;
		case 13:
			wbx.nucorg = "  " + n;
			break;
		default:
			wbx.nucorg = "FFF";
			break;
		}
		return wbx;
	}

	public Cell wbm2(int pi, int pj, char a, String n) {
		Cell wbx = m2(pi, pj, a, "" + Allookup.wb_lookup(n.charAt(0)).charAt(0) + Allookup.wb_lookup(n.charAt(1)).charAt(0));
		Cell wby = new Cell();
		for (int q = 0; q < Allookup.wb_lookup(n.charAt(0)).length(); q++) {
			for (int p = 0; p < Allookup.wb_lookup(n.charAt(1)).length(); p++) {
				wby = m2(pi, pj, a, "" + Allookup.wb_lookup(n.charAt(0)).charAt(q) + Allookup.wb_lookup(n.charAt(1)).charAt(p));
				if (mmMut(wby.score, wbx.score, this.wc_b_a)) {
					wbx = wby;
				}
			}
		}
		switch (wbx.casus) {
		case 21:
			wbx.nucorg = n + " ";
			break;
		case 22:
			wbx.nucorg = n.charAt(0) + " " + n.charAt(1);
			break;
		case 23:
			wbx.nucorg = " " + n;
			break;
		default:
			wbx.nucorg = "FFF";
			break;
		}
		return wbx;
	}

	public Cell wbm3(int pi, int pj, char a, String n) {
		Cell wbx = m3(pi, pj, a, "" + Allookup.wb_lookup(n.charAt(0)).charAt(0) + Allookup.wb_lookup(n.charAt(1)).charAt(0) + Allookup.wb_lookup(n.charAt(2)).charAt(0));
		Cell wby = new Cell();
		for (int q = 0; q < Allookup.wb_lookup(n.charAt(0)).length(); q++) {
			for (int p = 0; p < Allookup.wb_lookup(n.charAt(1)).length(); p++) {
				for (int r = 0; r < Allookup.wb_lookup(n.charAt(2)).length(); r++) {
					wby = m3(pi, pj, a, "" + Allookup.wb_lookup(n.charAt(0)).charAt(q) + Allookup.wb_lookup(n.charAt(1)).charAt(p) + Allookup.wb_lookup(n.charAt(2)).charAt(r));
					if (mmMut(wby.score, wbx.score, this.wc_b_a)) {
						wbx = wby;
					}
				}
			}
		}
		wbx.nucorg = n;
		return wbx;
	}

	public Cell wbm4(int pi, int pj, char a, String n) {
		Cell wbx = m4(pi, pj, a,
				"" + Allookup.wb_lookup(n.charAt(0)).charAt(0) + Allookup.wb_lookup(n.charAt(1)).charAt(0) + Allookup.wb_lookup(n.charAt(2)).charAt(0) + Allookup.wb_lookup(n.charAt(3)).charAt(0));
		Cell wby = new Cell();
		for (int q = 0; q < Allookup.wb_lookup(n.charAt(0)).length(); q++) {
			for (int p = 0; p < Allookup.wb_lookup(n.charAt(1)).length(); p++) {
				for (int r = 0; r < Allookup.wb_lookup(n.charAt(2)).length(); r++) {
					for (int s = 0; s < Allookup.wb_lookup(n.charAt(3)).length(); s++) {
						wby = m4(pi, pj, a, "" + Allookup.wb_lookup(n.charAt(0)).charAt(q) + Allookup.wb_lookup(n.charAt(1)).charAt(p) + Allookup.wb_lookup(n.charAt(2)).charAt(r)
								+ Allookup.wb_lookup(n.charAt(3)).charAt(s));
						if (mmMut(wby.score, wbx.score, this.wc_b_a)) {
							wbx = wby;
						}
					}
				}
			}
		}

		switch (wbx.casus) {
		case 41:
			wbx.nucorg = n.substring(0, 3);
			break;
		case 42:
			wbx.nucorg = n.substring(1, 4);
			break;
		case 43:
			wbx.nucorg = n.charAt(0) + n.substring(2, 4);
			break;
		case 44:
			wbx.nucorg = n.substring(0, 2) + n.charAt(3);
			break;
		default:
			wbx.nucorg = "FFF";
			break;
		}
		return wbx;
	}

	public Cell wbm5(int pi, int pj, char a, String n) {
		Cell wbx = m5(pi, pj, a,
				"" + Allookup.wb_lookup(n.charAt(0)).charAt(0) + Allookup.wb_lookup(n.charAt(1)).charAt(0) + Allookup.wb_lookup(n.charAt(2)).charAt(0) + Allookup.wb_lookup(n.charAt(3)).charAt(0)
						+ Allookup.wb_lookup(n.charAt(4)).charAt(0));
		Cell wby = new Cell();
		for (int q = 0; q < Allookup.wb_lookup(n.charAt(0)).length(); q++) {
			for (int p = 0; p < Allookup.wb_lookup(n.charAt(1)).length(); p++) {
				for (int r = 0; r < Allookup.wb_lookup(n.charAt(2)).length(); r++) {
					for (int s = 0; s < Allookup.wb_lookup(n.charAt(3)).length(); s++) {
						for (int t = 0; t < Allookup.wb_lookup(n.charAt(4)).length(); t++) {
							wby = m5(pi, pj, a, "" + Allookup.wb_lookup(n.charAt(0)).charAt(q) + Allookup.wb_lookup(n.charAt(1)).charAt(p) + Allookup.wb_lookup(n.charAt(2)).charAt(r)
									+ Allookup.wb_lookup(n.charAt(3)).charAt(s) + Allookup.wb_lookup(n.charAt(4)).charAt(t));
							if (mmMut(wby.score, wbx.score, this.wc_b_a)) {
								wbx = wby;
							}
						}
					}
				}
			}
		}
		switch (wbx.casus) {
		case 51:
			wbx.nucorg = n.substring(0, 3);
			break;
		case 52:
			wbx.nucorg = n.substring(0, 2) + n.charAt(3);
			break;
		case 53:
			wbx.nucorg = n.substring(0, 2) + n.charAt(4);
			break;
		case 54:
			wbx.nucorg = n.charAt(0) + n.substring(2, 4);
			break;
		case 55:
			wbx.nucorg = "" + n.charAt(0) + n.charAt(2) + n.charAt(4);
			break;
		case 56:
			wbx.nucorg = n.charAt(0) + n.substring(3);
			break;
		case 57:
			wbx.nucorg = n.substring(1, 4);
			break;
		case 58:
			wbx.nucorg = n.substring(1, 3) + n.charAt(4);
			break;
		case 59:
			wbx.nucorg = n.charAt(1) + n.substring(3);
			break;
		case 60:
			wbx.nucorg = n.substring(2);
			break;
		default:
			wbx.nucorg = "FFF";
			break;
		}

		return wbx;
	}

	public boolean wb_check(String wb) {
		Pattern r = Pattern.compile("[KSYMWRBDHVN]");
		Matcher m = r.matcher(wb);
		// System.out.println(m.find());
		return m.find();
	}

	public boolean mmMut(int wby, int wbx, boolean worstbest) {
		if (worstbest) {
			return (wbx <= -1000 || (wby < wbx && wby > -1000));
		} else {
			return (wby > wbx);
		}
	}

	public void align() { // old align without wild bases and stuff
		String p = sequence1;
		String n = sequence2.replaceAll("T", "U");

		s = new int[len1 + 1][len2 + 1];
		e = new int[len1 + 1][len2 + 1];
		h = new int[len1 + 1][len2 + 1];
		c_way = new Cell[len1 + 1][len2 + 1];
		Cell c1 = new Cell();
		Cell c2 = new Cell();
		Cell c3 = new Cell();
		Cell c4 = new Cell();
		Cell c5 = new Cell();
		for (int i = 0; i <= len1; i++) {
			c_way[i][0] = new Cell();
			s[i][0] = 0;
			e[i][0] = 0;
			h[i][0] = 0;
		}
		for (int j = 0; j <= len2; j++) {
			c_way[0][j] = new Cell();
			s[0][j] = 0;
			e[0][j] = 0;
			h[0][j] = 0;
		}

		for (int i = 1; i < len1; i++) {
			for (int j = 5; j < len2; j++) {
				e[i][j] = Math.max(s[i - 1][j] + gop, e[i - 1][j] + gep);
				h[i][j] = Math.max(s[i][j - 3] + gop, h[i][j - 3] + gep);

				c1 = m1(i - 1, j - 1, p.charAt(i), n.charAt(j));
				c2 = m2(i - 1, j - 2, p.charAt(i), "" + n.charAt(j - 1) + n.charAt(j));
				c3 = m3(i - 1, j - 3, p.charAt(i), "" + n.charAt(j - 2) + n.charAt(j - 1) + n.charAt(j));
				c4 = m4(i - 1, j - 4, p.charAt(i), "" + n.charAt(j - 3) + n.charAt(j - 2) + n.charAt(j - 1) + n.charAt(j));
				c5 = m5(i - 1, j - 5, p.charAt(i), "" + n.charAt(j - 4) + n.charAt(j - 3) + n.charAt(j - 2) + n.charAt(j - 1) + n.charAt(j));

				s[i][j] = Math.max(Util.max(s[i - 1][j - 1] + c1.score, s[i - 1][j - 2] + c2.score, s[i - 1][j - 3] + c3.score, s[i - 1][j - 4] + c4.score),
						Util.max(s[i - 1][j - 5] + c5.score, e[i][j], h[i][j], 0));
				if (s[i][j] == s[i - 1][j - 1] + c1.score) {
					c_way[i][j] = c1;
				}
				if (s[i][j] == s[i - 1][j - 2] + c2.score) {
					c_way[i][j] = c2;
				}
				if (s[i][j] == s[i - 1][j - 3] + c3.score) {
					c_way[i][j] = c3;
				}
				if (s[i][j] == s[i - 1][j - 4] + c4.score) {
					c_way[i][j] = c4;
				}
				if (s[i][j] == s[i - 1][j - 5] + c5.score) {
					c_way[i][j] = c5;
				}
				if (s[i][j] == e[i][j]) {
					c_way[i][j] = new Cell();
					c_way[i][j].casus = 7;
					c_way[i][j].p_row = i - 1;
					c_way[i][j].p_col = j;
					c_way[i][j].acid = p.charAt(i);
					c_way[i][j].nuc = "---";
					c_way[i][j].triplet = "---";
				}
				if (s[i][j] == h[i][j]) {
					c_way[i][j] = new Cell();
					c_way[i][j].casus = 8;
					c_way[i][j].p_row = i;
					c_way[i][j].p_col = j - 3;
					c_way[i][j].acid = '_';
					c_way[i][j].nuc = "" + n.charAt(j - 3) + n.charAt(j - 2) + n.charAt(j - 1);
					c_way[i][j].triplet = "" + n.charAt(j - 3) + n.charAt(j - 2) + n.charAt(j - 1);
				}
				if (s[i][j] == 0) {
					c_way[i][j] = new Cell();
					c_way[i][j].casus = -42;
				}

				if (s[i][j] > bestval) {
					bestval = s[i][j];
					bestrow = i; // row;
					bestcol = j;
				}
			}
		}
		System.out.println("alignment score:" + bestval);
	}

	public Cell m1(int pi, int pj, char a, char n) {
		// XAA, AXA, AAX, usw....also 2 inserts
		Cell x = new Cell();

		int max = -100000;
		int perm1 = -1, perm2 = -1;
		int f = -1;
		char[] l = { 'A', 'C', 'G', 'U' };
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				if (max < matrixCost(Allookup.lookup(a), Allookup.lookup(Allookup.trilookup("" + n + l[i] + l[j])), sc_matrix)) {
					max = matrixCost(Allookup.lookup(a), Allookup.lookup(Allookup.trilookup("" + n + l[i] + l[j])), sc_matrix);
					perm1 = i;
					perm2 = j;
					f = 1;
				}
				if (max < matrixCost(Allookup.lookup(a), Allookup.lookup(Allookup.trilookup("" + l[i] + n + l[j])), sc_matrix)) {
					max = matrixCost(Allookup.lookup(a), Allookup.lookup(Allookup.trilookup("" + l[i] + n + l[j])), sc_matrix);
					perm1 = i;
					perm2 = j;
					f = 2;
				}
				if (max < matrixCost(Allookup.lookup(a), Allookup.lookup(Allookup.trilookup("" + l[i] + l[j] + n)), sc_matrix)) {
					max = matrixCost(Allookup.lookup(a), Allookup.lookup(Allookup.trilookup("" + l[i] + l[j] + n)), sc_matrix);
					perm1 = i;
					perm2 = j;
					f = 3;
				}
			}
		}
		x.val = perm1 * 100 + perm2 + f;

		if (f == -1) {
			x.triplet = "F1F";
		}
		if (f == 1) {
			x.triplet = "" + n + l[perm1] + l[perm2];
			x.nuc = n + "  ";
			x.casus = 11;
		}
		if (f == 2) {
			x.triplet = "" + l[perm1] + n + l[perm2];
			x.nuc = " " + n + " ";
			x.casus = 12;
		}
		if (f == 3) {
			x.triplet = "" + l[perm1] + l[perm2] + n;
			x.nuc = "  " + n;
			x.casus = 13;
		}
		x.acid = a;
		x.score = max + dfsp;
		x.p_row = pi;
		x.p_col = pj;
		x.nucorg = x.triplet;
		return x;
	}

	public Cell m2(int pi, int pj, char a, String n) {
		// AAX AXA XAA , also 1 insert
		Cell x = new Cell();
		int max = -1000000;
		int perm1 = -1;
		int f = -1;
		char[] l = { 'A', 'C', 'G', 'U' };
		for (int i = 0; i < 4; i++) {
			if (max < matrixCost(Allookup.lookup(a), Allookup.lookup(Allookup.trilookup(n + l[i])), sc_matrix)) {
				max = matrixCost(Allookup.lookup(a), Allookup.lookup(Allookup.trilookup(n + l[i])), sc_matrix);
				perm1 = i;
				f = 1;
			}
			if (max < matrixCost(Allookup.lookup(a), Allookup.lookup(Allookup.trilookup("" + n.charAt(0) + l[i] + n.charAt(1))), sc_matrix)) {
				max = matrixCost(Allookup.lookup(a), Allookup.lookup(Allookup.trilookup("" + n.charAt(0) + l[i] + n.charAt(1))), sc_matrix);
				perm1 = i;
				f = 2;
			}
			if (max < matrixCost(Allookup.lookup(a), Allookup.lookup(Allookup.trilookup(l[i] + n)), sc_matrix)) {
				max = matrixCost(Allookup.lookup(a), Allookup.lookup(Allookup.trilookup(l[i] + n)), sc_matrix);
				perm1 = i;
				f = 3;
			}
		}
		x.val = perm1 * 10 + f;
		if (f == -1) {
			x.triplet = "F_F";
		}
		if (f == 1) {
			x.triplet = n + l[perm1];
			x.nuc = n + " ";
			x.casus = 21;
		}
		if (f == 2) {
			x.triplet = "" + n.charAt(0) + l[perm1] + n.charAt(1);
			x.nuc = n.charAt(0) + " " + n.charAt(1);
			x.casus = 22;
		}
		if (f == 3) {
			x.triplet = l[perm1] + n;
			x.nuc = " " + n;
			x.casus = 23;
		}

		x.acid = a;
		x.score = max + sfsp;
		x.p_row = pi;
		x.p_col = pj;
		// x.nucorg = x.nuc;
		return x;
	}

	public Cell m3(int pi, int pj, char a, String n) {
		// AAA
		Cell x = new Cell();
		x.acid = a;
		x.triplet = n;
		x.nuc = n;
		x.nucorg = x.nuc;
		x.val = s[pi][pj];
		x.p_row = pi;
		x.p_col = pj;
		x.casus = 3;
		x.score = matrixCost(Allookup.lookup(a), Allookup.lookup(Allookup.trilookup(n)), sc_matrix);

		return x;
	}

	public Cell m4(int pi, int pj, char a, String n) {
		// AAAA und eins raus (1 del)
		Cell x = new Cell();
		int v1 = matrixCost(Allookup.lookup(a), Allookup.lookup(Allookup.trilookup(n.substring(0, 3))), sc_matrix); // 012
		int v2 = matrixCost(Allookup.lookup(a), Allookup.lookup(Allookup.trilookup(n.substring(1, 4))), sc_matrix); // 123
		int v3 = matrixCost(Allookup.lookup(a), Allookup.lookup(Allookup.trilookup(n.charAt(0) + n.substring(2, 4))), sc_matrix); // 023
		int v4 = matrixCost(Allookup.lookup(a), Allookup.lookup(Allookup.trilookup(n.substring(0, 2) + n.charAt(3))), sc_matrix); // 013

		if (v1 >= v2 && v1 >= v3 && v1 >= v4) {
			x.score = v1;
			x.triplet = n.substring(0, 3);
			// anpassung org-Sequenz für ausgabe via nucorg
			x.nucorg = n.substring(0, 3);
			x.nuc = n;
			x.casus = 41;
		}
		if (v2 >= v1 && v2 >= v3 && v2 >= v4) {
			x.score = v2;
			x.triplet = n.substring(1, 4);
			// anpassung org-Sequenz für ausgabe via nucorg
			x.nucorg = n.substring(1, 4);
			x.nuc = n;
			x.casus = 42;
		}
		if (v3 >= v1 && v3 >= v2 && v3 >= v4) {
			x.score = v3;
			x.triplet = n.charAt(0) + n.substring(2, 4);
			// anpassung org-Sequenz für ausgabe via nucorg
			x.nucorg = n.charAt(0) + n.substring(2, 4);
			x.nuc = n;
			x.casus = 43;
		}
		if (v4 >= v1 && v4 >= v2 && v4 >= v3) {
			x.score = v4;
			x.triplet = n.substring(0, 2) + n.charAt(3);
			// anpassung org-Sequenz für ausgabe via nucorg
			x.nucorg = n.substring(0, 2) + n.charAt(3);
			x.nuc = n;
			x.casus = 44;
		}

		x.score = x.score + sfsp;
		x.acid = a;
		x.val = s[pi][pj];
		x.p_row = pi;
		x.p_col = pj;
		// x.nucorg = x.nuc;
		return x;
	}

	public Cell m5(int pi, int pj, char a, String n) {
		// AAAAA, also 2 del
		Cell x = new Cell();
		int[] vals = new int[10];
		vals[0] = matrixCost(Allookup.lookup(a), Allookup.lookup(Allookup.trilookup(n.substring(0, 3))), sc_matrix);
		vals[1] = matrixCost(Allookup.lookup(a), Allookup.lookup(Allookup.trilookup(n.substring(0, 2) + n.charAt(3))), sc_matrix);
		vals[2] = matrixCost(Allookup.lookup(a), Allookup.lookup(Allookup.trilookup(n.substring(0, 2) + n.charAt(4))), sc_matrix);
		vals[3] = matrixCost(Allookup.lookup(a), Allookup.lookup(Allookup.trilookup(n.charAt(0) + n.substring(2, 4))), sc_matrix);
		vals[4] = matrixCost(Allookup.lookup(a), Allookup.lookup(Allookup.trilookup("" + n.charAt(0) + n.charAt(2) + n.charAt(4))), sc_matrix);
		vals[5] = matrixCost(Allookup.lookup(a), Allookup.lookup(Allookup.trilookup(n.charAt(0) + n.substring(3))), sc_matrix);
		vals[6] = matrixCost(Allookup.lookup(a), Allookup.lookup(Allookup.trilookup(n.substring(1, 4))), sc_matrix);
		vals[7] = matrixCost(Allookup.lookup(a), Allookup.lookup(Allookup.trilookup(n.substring(1, 3) + n.charAt(4))), sc_matrix);
		vals[8] = matrixCost(Allookup.lookup(a), Allookup.lookup(Allookup.trilookup(n.charAt(1) + n.substring(3))), sc_matrix);
		vals[9] = matrixCost(Allookup.lookup(a), Allookup.lookup(Allookup.trilookup(n.substring(2))), sc_matrix);

		x.score = Util.max(Util.max(vals[0], vals[1], vals[2], vals[3]), Util.max(vals[4], vals[5], vals[6], vals[7]), vals[8], vals[9]);
		// System.out.println(x.score);
		for (int i = 0; i < 10; i++) {
			if (x.score == vals[i]) {
				switch (i) {
				case 0: {
					x.triplet = n.substring(0, 3);
					x.nucorg = n.substring(0, 3);
					x.nuc = n;
					x.casus = 51;
					break;
				}
				case 1: {
					x.triplet = n.substring(0, 2) + n.charAt(3);
					x.nucorg = n.substring(0, 2) + n.charAt(3);
					x.nuc = n;
					x.casus = 52;
					break;
				}
				case 2: {
					x.triplet = n.substring(0, 2) + n.charAt(4);
					x.nucorg = n.substring(0, 2) + n.charAt(4);
					x.nuc = n;
					x.casus = 53;
					break;
				}
				case 3: {
					x.triplet = n.charAt(0) + n.substring(2, 4);
					x.nucorg = n.charAt(0) + n.substring(2, 4);
					x.nuc = n;
					x.casus = 54;
					break;
				}
				case 4: {
					x.triplet = "" + n.charAt(0) + n.charAt(2) + n.charAt(4);
					x.nucorg = "" + n.charAt(0) + n.charAt(2) + n.charAt(4);
					x.nuc = n;
					x.casus = 55;
					break;
				}
				case 5: {
					x.triplet = n.charAt(0) + n.substring(3);
					x.nucorg = n.charAt(0) + n.substring(3);
					x.nuc = n;
					x.casus = 56;
					break;
				}
				case 6: {
					x.triplet = n.substring(1, 4);
					x.nucorg = n.substring(1, 4);
					x.nuc = n;
					x.casus = 57;
					break;
				}
				case 7: {
					x.triplet = n.substring(1, 3) + n.charAt(4);
					x.nucorg = n.substring(1, 3) + n.charAt(4);
					x.nuc = n;
					x.casus = 58;
					break;
				}
				case 8: {
					x.triplet = n.charAt(1) + n.substring(3);
					x.nucorg = n.charAt(1) + n.substring(3);
					x.nuc = n;
					x.casus = 59;
					break;
				}
				case 9: {
					x.triplet = n.substring(2);
					x.nucorg = n.substring(2);
					x.nuc = n;
					x.casus = 60;
					break;
				}
				default: {
					x.triplet = "F2F";
					x.nuc = n;
					break;
				}
				}
				break;
			}
		}
		x.score = x.score + dfsp;
		x.acid = a;
		x.val = s[pi][pj];
		x.p_row = pi;
		x.p_col = pj;
		// x.nucorg = x.nuc;
		return x;
	}

	public int matrixCost(int i, int j, int[][] sc_matrix) {
		return sc_matrix[i][j];
	}

	public void ausgabe(int[][] s) {
		for (int i = 0; i < s.length; i++) {
			for (int j = 0; j < s[0].length; j++) {
				System.out.print(" " + s[i][j] + " ,");
			}
			System.out.println("");
		}
	}

	public void ausgCell(Cell[][] c) {
		for (int i = 0; i < c.length; i++) {
			for (int j = 0; j < c[0].length; j++) {
				if (c[i][j] == null) {
					System.out.print("-+- ,");
				} else {
					System.out.print(" " + c[i][j].p_col + "+" + c[i][j].p_row + " ,");
				}
			}
			System.out.println("");
		}
	}

	public void showWay() {
		int i = bestrow;
		int j = bestcol;
		System.out.print("von" + bestrow);
		System.out.println("|" + bestcol);
		while (s[i][j] > 0) {
			System.out.print(c_way[i][j].acid + "--> : ");
			System.out.print(c_way[i][j].nuc + "  go to: ");
			System.out.print(c_way[i][j].p_row + "|" + c_way[i][j].p_col);
			System.out.println(" with score: " + s[i][j]);
			int m = c_way[i][j].p_row;
			j = c_way[i][j].p_col;
			i = m;
			// System.out.println(i +"->"+ j);
			if (i < 0 || j < 0) {
				break;
			}
		}
	}

	public void setAlignment() { // Zwischenmarkierungen setzten + Sequenzen
		int i = bestrow;
		int j = bestcol;
		StringBuffer p = new StringBuffer("");
		StringBuffer a = new StringBuffer("");
		StringBuffer n = new StringBuffer("");
		StringBuffer o = new StringBuffer("");
		StringBuffer nfix = new StringBuffer("");

		while (s[i][j] > 0) {
			switch (c_way[i][j].casus) {
			case 11: {
				a.insert(0, "|ii");
				o.insert(0, " ii");
				break;
			}
			case 12: {
				a.insert(0, "i|i");
				o.insert(0, "i i");
				break;
			}
			case 13: {
				a.insert(0, "ii|");
				o.insert(0, "ii ");
				break;
			}
			case 21: {
				a.insert(0, "||i");
				o.insert(0, "  i");
				break;
			}
			case 22: {
				a.insert(0, "|i|");
				o.insert(0, " i ");
				break;
			}
			case 23: {
				a.insert(0, "i||");
				o.insert(0, "i  ");
				break;
			}
			case 3: {
				if (c_way[i][j].acid == Allookup.trilookup(c_way[i][j].nuc)) {
					a.insert(0, "|||");
					o.insert(0, "   ");
					mut.idens++;
				} else {
					if (c_way[i][j].score >= 0) {
						a.insert(0, "+++");
						o.insert(0, "   ");
						mut.posis++;
					} else {
						a.insert(0, "***");
						o.insert(0, "RRR");
					}
				}
				break;
			}
			case 41: {
				a.insert(0, "|-4");
				o.insert(0, "   d");
				break;
			}
			case 42: {
				a.insert(0, "|-1");
				o.insert(0, "d   ");
				break;
			}
			case 43: {
				a.insert(0, "|-2");
				o.insert(0, " d  ");
				break;
			}
			case 44: {
				a.insert(0, "|-3");
				o.insert(0, "  d ");
				break;
			}
			case 51: {
				a.insert(0, "-45");
				o.insert(0, "   dd");
				break;
			}
			case 52: {
				a.insert(0, "-35");
				o.insert(0, "  d d");
				break;
			}
			case 53: {
				a.insert(0, "-34");
				o.insert(0, "  dd ");
				break;
			}
			case 54: {
				a.insert(0, "-25");
				o.insert(0, " d  d");
				break;
			}
			case 55: {
				a.insert(0, "-24");
				o.insert(0, " d d ");
				break;
			}
			case 56: {
				a.insert(0, "-23");
				o.insert(0, "  dd ");
				break;
			}
			case 57: {
				a.insert(0, "-15");
				o.insert(0, "d   d");
				break;
			}
			case 58: {
				a.insert(0, "-14");
				o.insert(0, "d  d ");
				break;
			}
			case 59: {
				a.insert(0, "-13");
				o.insert(0, "d d  ");
				break;
			}
			case 60: {
				a.insert(0, "-12");
				o.insert(0, "dd   ");
				break;
			}
			case 7: {
				a.insert(0, "DDD");
				o.insert(0, "DDD");
				break;
			}
			case 8: {
				a.insert(0, "III");
				o.insert(0, "III");
				break;
			}
			case 9: {
				a.insert(0, "STP");
				o.insert(0, "STP");
				break;
			}
			default: {
				a.insert(0, "???");
				o.insert(0, "???");
				break;
			}
			}
			p.insert(0, Allookup.p2pro(c_way[i][j].acid));

			if (c_way[i][j].nucorg.length() > 0) {
				n.insert(0, c_way[i][j].nucorg);
			} else {
				n.insert(0, c_way[i][j].nuc);
			}

			nfix.insert(0, c_way[i][j].triplet);
			int m = c_way[i][j].p_row;
			j = c_way[i][j].p_col;
			i = m;
		}

		seq1_al = p.toString();
		alig_al = a.toString();
		seq2o_al = n.toString();
		alig_o = o.toString();
		seq2_al = nfix.toString();
	}

	public void showAlignment() {
		setAlignment();
		System.out.println();
		int r1 = startrow;
		int r2 = startcol;
		int res = this.res;
		for (int i = 0; i < seq1_al.length(); i += res) {

			if (i + res > seq1_al.length()) {
				res = seq1_al.length() - i;
			}
			System.out.printf("%6s", r1);
			System.out.println(" " + seq1_al.substring(i, i + res));
			System.out.printf("%7s", "");
			System.out.println(alig_al.substring(i, i + res));
			System.out.printf("%6s", r2);
			System.out.println(" " + seq2_al.substring(i, i + res));

			// MARKER
			// System.out.printf("%7s", "");
			// System.out.println(alig_o.substring(i, i + res));

			// orginal Sequence....[DEBUG]
			// System.out.printf("%7s", "");
			// System.out.println(seq2o_al.substring(i, i + res));
			// SEQUENCE!!....[DEBUG]
			// System.out.printf("%4s", "");
			// System.out.println(sequence2);

			System.out.println("");
			r1 += res / 3;
			r2 += res;
		}
		/*
		 * System.out.println(seq1_al); System.out.println(alig_al);
		 * System.out.println(seq2o_al); System.out.println(seq2_al);
		 */
	}

	public void setSize() {
		int i = bestrow;
		int j = bestcol;
		int l = 0;
		while (s[i][j] > 0) {
			l++;
			int m = c_way[i][j].p_row;
			j = c_way[i][j].p_col;
			i = m;
			if (i < 0 || j < 0) {
				break;
			}
		}
		startrow = i - 1; // -2;
		startcol = j - 1; // -3;
		len_al = l;
	}

	public void showSize() {
		System.out.print("Alignment: from " + startrow + "[p]|" + startcol + "[n]  ");
		System.out.println("to " + (bestrow - 2) + "[p]|" + (bestcol - 3) + "[n] ");

	}

	public void setMut() { // ANPASSUNG AUCH FÜR JSON-FILE! ! ! !
		mut = new Mutator();
		int i = bestrow;
		int j = bestcol;
		int mc1 = 0;
		int mc1p = 0;
		int mc2 = 0;
		int mc3 = 0;
		int mc4_1 = 0;
		int mc4_2 = 0;
		StringBuffer i_pos = new StringBuffer("");
		StringBuffer ix_pos = new StringBuffer("");
		StringBuffer d_pos = new StringBuffer("");
		StringBuffer dx_pos = new StringBuffer("");
		StringBuffer r_pos = new StringBuffer("");
		StringBuffer rx_pos = new StringBuffer("");
		StringBuffer px_pos = new StringBuffer("");
		StringBuffer sf_pos = new StringBuffer("");
		StringBuffer df_pos = new StringBuffer("");
		StringBuffer sfx_pos = new StringBuffer("");
		StringBuffer dfx_pos = new StringBuffer("");
		// ArrayList<Tuple> in_pos = new ArrayList<Tuple>();

		while (s[i][j] > 0) { // nuc durch nucorg btw. triplet durch nucorg
								// ersetzt
			switch (c_way[i][j].casus) {
			case 11: {
				mc4_2++;
				break;
			}
			case 12: {
				mc4_2++;
				break;
			}
			case 13: {
				mc4_2++;
				break;
			}
			case 21: {
				mc4_1++;
				break;
			}
			case 22: {
				mc4_1++;
				break;
			}
			case 23: {
				mc4_1++;
				break;
			}
			case 3: { // REPLACE or POS
				String wbp = "";
				if(wb_check(c_way[i][j].nucorg)){
					wbp = "["+c_way[i][j].nuc + "]->"+ Allookup.trilookup(c_way[i][j].nuc);
				} else{
					wbp = "->"+ Allookup.trilookup(c_way[i][j].nucorg);
				}
				if (c_way[i][j].score < 0) {
					mc1++;
					r_pos.insert(0, (c_way[i][j].p_col - 1) + "  ");
		//			rx_pos.insert(0, (c_way[i][j].p_row - 1) + "(" + c_way[i][j].acid + ":" + Allookup.trilookup(c_way[i][j].nucorg) + ")  ");
					rx_pos.insert(0, (c_way[i][j].p_row - 1) + "(" + c_way[i][j].acid + ":" + c_way[i][j].nucorg + wbp +")  ");
					mut.re_pos.add(new Tuple(c_way[i][j].p_row - 1, 1));
				} else {
					if (c_way[i][j].acid != Allookup.trilookup(c_way[i][j].nuc)) {
						mc1p++;
		//				px_pos.insert(0, (c_way[i][j].p_row - 1) + "(" + c_way[i][j].acid + ":" + Allookup.trilookup(c_way[i][j].nucorg) + ")  ");
						px_pos.insert(0, (c_way[i][j].p_row - 1) + "(" + c_way[i][j].acid + ":" + c_way[i][j].nucorg + wbp + ")  ");
						mut.po_pos.add(new Tuple(c_way[i][j].p_row - 1, 1));
					}
				}
				break;
			}
			case 41: {
				mc4_1++;
				break;
			}
			case 42: {
				mc4_1++;
				break;
			}
			case 43: {
				mc4_1++;
				break;
			}
			case 44: {
				mc4_1++;
				break;
			}
			case 51: {
				mc4_2++;
				break;
			}
			case 52: {
				mc4_2++;
				break;
			}
			case 53: {
				mc4_2++;
				break;
			}
			case 54: {
				mc4_2++;
				break;
			}
			case 55: {
				mc4_2++;
				break;
			}
			case 56: {
				mc4_2++;
				break;
			}
			case 57: {
				mc4_2++;
				break;
			}
			case 58: {
				mc4_2++;
				break;
			}
			case 59: {
				mc4_2++;
				break;
			}
			case 60: {
				mc4_2++;
				break;
			}
			case 7: { // DELETION
				mc3++;
				d_pos.insert(0, (c_way[i][j].p_col - 1) + "  ");
				dx_pos.insert(0, (c_way[i][j].p_row - 1) + "(-:" + c_way[i][j].acid + ")  ");
				int runner = 0;
				boolean found = false;
				while (runner < mut.de_pos.size()) {
					if (mut.de_pos.get(runner).index == (c_way[i][j].p_row - 1)) {
						mut.de_pos.set(runner, new Tuple(mut.de_pos.get(runner).index, (mut.de_pos.get(runner).count) + 1));
						found = true;
						break;
					}
					runner++;
				}
				if (!found) {
					mut.de_pos.add(new Tuple((c_way[i][j].p_row - 1), 1));
				}
				break;
			}
			case 8: { // INSERTION
				mc2++;
				i_pos.insert(0, (c_way[i][j].p_col - 1) + "  ");
				ix_pos.insert(0, (c_way[i][j].p_row - 1) + "(" + Allookup.trilookup(c_way[i][j].nucorg) + ":-)  ");
				int runner = 0;
				boolean found = false;
				while (runner < mut.in_pos.size()) {
					if (mut.in_pos.get(runner).index == (c_way[i][j].p_row - 1)) {
						mut.in_pos.set(runner, new Tuple(mut.in_pos.get(runner).index, (mut.in_pos.get(runner).count) + 1));
						found = true;
						break;
					}
					runner++;
				}
				if (!found) {
					mut.in_pos.add(new Tuple((c_way[i][j].p_row - 1), 1));
				}
				break;
			}
			default:
				break;
			}

			switch (c_way[i][j].casus) {
			case 11:
			case 12:
			case 13:
			case 51:
			case 52:
			case 53:
			case 54:
			case 55:
			case 56:
			case 57:
			case 58:
			case 59:
			case 60: { // DOUBLE FRAME SHIFT
				df_pos.insert(0, (c_way[i][j].p_col - 1) + "  ");
				// dfx_pos.insert(0, (c_way[i][j].p_row - 1) + "(" +
				// c_way[i][j].acid + ":" +
				// Allookup.trilookup(c_way[i][j].triplet) + ")  ");
				dfx_pos.insert(0, (c_way[i][j].p_row - 1) + "(DFS:" + Allookup.trilookup(c_way[i][j].nucorg) + ")  ");
				mut.do_pos.add(new Tuple((c_way[i][j].p_row - 1), 1));
				break;
			}
			case 21:
			case 22:
			case 23:
			case 41:
			case 42:
			case 43:
			case 44: { // SINGLE FRAME SHIFT
				sf_pos.insert(0, (c_way[i][j].p_col - 1) + "  ");
				// sfx_pos.insert(0, (c_way[i][j].p_row - 1) + "(" +
				// c_way[i][j].acid + ":" +
				// Allookup.trilookup(c_way[i][j].triplet) + ")  ");
				sfx_pos.insert(0, (c_way[i][j].p_row - 1) + "(SFS:" + Allookup.trilookup(c_way[i][j].nucorg) + ")  ");
				mut.si_pos.add(new Tuple((c_way[i][j].p_row - 1), 1));
				break;
			}
			default:
				break;
			}

			int m = c_way[i][j].p_row;
			j = c_way[i][j].p_col;
			i = m;
			if (i < 0 || j < 0) {
				break;
			}
		}
		mut.inserts = mc2;
		mut.deletions = mc3;
		mut.positives = mc1p;
		mut.replaces = mc1;
		mut.singlefs = mc4_1;
		mut.doublefs = mc4_2;
		mut.frameshifts = mc4_1 + mc4_2;

		// mut.i_pos = i_pos.toString();
		mut.i_pos = ix_pos.toString();
		// try for graph
		// System.out.println(mut.i_pos);
		// try for graph

		// mut.d_pos = d_pos.toString();
		mut.d_pos = dx_pos.toString();
		// mut.r_pos = r_pos.toString();
		mut.p_pos = px_pos.toString();
		mut.r_pos = rx_pos.toString();
		mut.sf_pos = sf_pos.toString();
		mut.sf_pos = sfx_pos.toString();
		mut.df_pos = df_pos.toString();
		mut.df_pos = dfx_pos.toString();
	}

	public void showMut() {
		System.out.println("Mutations:");
		System.out.print("# of INS: " + mut.inserts);
		System.out.println(" [positions of INS: " + mut.i_pos + "]");
		System.out.print("# of DEL: " + mut.deletions);
		System.out.println(" [positions of DEL: " + mut.d_pos + "]");
		System.out.print("# of REP: " + mut.replaces);
		System.out.println(" [positions of REP: " + mut.r_pos + "]");
		System.out.print("# of FSM: " + mut.frameshifts + " ( sfs: " + mut.singlefs + "  and dfs: " + mut.doublefs + ") ");
		System.out.println(" [positions of SFS: " + mut.sf_pos + "] and [positions of DFS: " + mut.df_pos + "]");
	}

	public void showTuple() { // TESTING ONLY!
		System.out.println("Tuples:");
		System.out.println("INS: " + mut.in_pos);
		System.out.println("DEL: " + mut.de_pos);
		System.out.println("REP: " + mut.re_pos);
		System.out.println("POS: " + mut.po_pos);
		System.out.println("SFS: " + mut.si_pos);
		System.out.println("DFS: " + mut.do_pos);
		System.out.println("STP: " + mut.st_pos);
	}

	public void setStop() {
		int i = bestrow;
		int j = bestcol;
		StringBuffer stp_pos = new StringBuffer("");
		StringBuffer stx_pos = new StringBuffer("");

		while (s[i][j] > 0) {
			switch (c_way[i][j].triplet) {
			case "UAA":
			case "UGA":
			case "UAG":
			case "TAA":
			case "TGA":
			case "TAG": {
				stp_pos.insert(0, "[" + ((c_way[i][j].p_col - 1) + " to " + (c_way[i][j].p_col - 1 + 3) + "]  "));
				stx_pos.insert(0, (c_way[i][j].p_row - 1) + "(Stop:" + Allookup.trilookup(c_way[i][j].triplet) + ")  ");
				mut.st_pos.add(new Tuple(c_way[i][j].p_row - 1, 1));
			}

			default:
			}
			int m = c_way[i][j].p_row;
			j = c_way[i][j].p_col;
			i = m;
			if (i < 0 || j < 0) {
				break;
			}
		}
		mut.stp_pos = stx_pos.toString();
	}

	public void setWB() {
		int[] wb = new int[11];
		for (int i = 0; i < sequence2.length(); i++) {
			switch (sequence2.charAt(i)) {
			case 'K':
				wb[0]++;
				break;
			case 'S':
				wb[1]++;
				break;
			case 'Y':
				wb[2]++;
				break;
			case 'M':
				wb[3]++;
				break;
			case 'W':
				wb[4]++;
				break;
			case 'R':
				wb[5]++;
				break;
			case 'B':
				wb[6]++;
				break;
			case 'D':
				wb[7]++;
				break;
			case 'H':
				wb[8]++;
				break;
			case 'V':
				wb[9]++;
				break;
			case 'N':
				wb[10]++;
				break;
			default:
				break;
			}
		}
		wb_counter = wb;
	}

	public void showStop() {
		if (mut.stp_pos.length() > 0) {
			System.out.println("Warning!");
			System.out.println("detected STOP-codon in the alignment at " + mut.stp_pos);
		}
	}

	public void showNames() {
		System.out.println("Protein sequence: " + name1);
		System.out.println("Nucleotide sequence: " + name2);
	}

	public void showScore() {
		System.out.println("Alignment score: " + bestval);
	}

	public void showGaps() {
		NumberFormat formatter = new DecimalFormat("00.00");
		System.out.println("# of gaps: " + (mut.inserts + mut.deletions) + "[# of INS:" + mut.inserts + " and # of DELS: " + mut.deletions + "]" + "("
				+ formatter.format(((mut.inserts + mut.deletions) * 100.0 / (bestrow - 2 - startrow))) + "%)");
	}

	public void showIdens() {
		NumberFormat formatter = new DecimalFormat("00.00");
		System.out.println("# of identical matches: " + mut.idens + "(" + formatter.format((mut.idens * 100.0 / (bestrow - 2 - startrow))) + "%)");
	}

	public void showPosis() {
		NumberFormat formatter = new DecimalFormat("00.00");
		System.out.println("# of positive matches: " + mut.posis + "(" + formatter.format((mut.posis * 100.0 / (bestrow - 2 - startrow))) + "%)");
	}

	public void showSLength() {
		System.out.println("Length of the protein sequence:" + (len1 - 3));
		System.out.println("Length of the nucleotide sequence:" + (len2 - 3));
	}

	public void showALength() {
		System.out.println("Length of the alignment:" + (bestrow - 2 - startrow) + "[" + (len_al + mut.deletions - mut.inserts) + "]");
	}

	public void showParameters() {
		System.out.print("Parameters: ");
		System.out.print("gop: " + gop + ", ");
		System.out.print("gep: " + gep + ", ");
		System.out.print("sfsp: " + sfsp + ", ");
		System.out.print("dfsp: " + dfsp + ", ");
		if (matrix.contains("/")) {
			System.out.print("matrix: " + matrix.substring(matrix.lastIndexOf('/') + 1) + ", ");
		} else {
			System.out.print("matrix: " + matrix + ", ");
		}
		System.out.println("special: " + special);
	}

	public void showWBs() {
		setWB();
		System.out.print("# of wild bases: ");
		System.out
				.print(wb_counter[0] + wb_counter[1] + wb_counter[2] + wb_counter[3] + wb_counter[4] + wb_counter[5] + wb_counter[6] + wb_counter[7] + wb_counter[8] + wb_counter[9] + wb_counter[10]);
		System.out.print(" ( K=" + wb_counter[0] + ", ");
		System.out.print("S=" + wb_counter[1] + ", ");
		System.out.print("Y=" + wb_counter[2] + ", ");
		System.out.print("M=" + wb_counter[3] + ", ");
		System.out.print("W=" + wb_counter[4] + ", ");
		System.out.print("R=" + wb_counter[5] + ", ");
		System.out.print("B=" + wb_counter[6] + ", ");
		System.out.print("D=" + wb_counter[7] + ", ");
		System.out.print("H=" + wb_counter[8] + ", ");
		System.out.print("V=" + wb_counter[9] + ", ");
		System.out.print("N=" + wb_counter[10] + " )");
		System.out.println();
	}

	public static String[] readFile(String file) throws IOException {
		String[] s_n = new String[2];
		String inhalt = "";
		String zeile = "";
		StringBuffer sb = new StringBuffer();
		BufferedReader in = new BufferedReader(new FileReader(file));
		zeile = in.readLine();
		if (zeile.startsWith(">")) {
			s_n[0] = zeile.substring(1);
		} else {
			s_n[0] = "unnamed sequence";
			inhalt = inhalt + zeile;
		}
		while ((zeile = in.readLine()) != null) {
			inhalt = inhalt + zeile.replaceAll("[\r\n]+", "").trim();
		}
		in.close();
		s_n[1] = sb.append(inhalt).toString();
		// System.out.println(s_n[1]);
		return s_n;
	}

	public static MString jsonp(String source) {
		int offset = 0;
		MString parts = new MString();
		String[] tokens;
		if (source.length() > 1) {
			tokens = source.split("[ ]+");
			// if (tokens.length >= 1) {
			for (int z = 0; z < tokens.length - 1; z++) {
				parts.s1 = parts.s1 + "{\"position\":" + tokens[z].substring(0, tokens[z].indexOf("(")) + ",\"subject\":\"" + tokens[z].substring(tokens[z].indexOf("(") + 1, tokens[z].indexOf(":"))
						+ "\",\"query\":\"" + tokens[z].substring(tokens[z].indexOf(":") + 1, tokens[z].indexOf(")")) + "\"},";
				parts.s2 = parts.s2 + "{\"begin\":" + tokens[z].substring(0, tokens[z].indexOf("(")) + ",\"end\":" + (Integer.parseInt(tokens[z].substring(0, tokens[z].indexOf("("))) + offset);
				parts.s2 = parts.s2 + "},";
			}
			parts.s1 = parts.s1 + "{";
			parts.s1 = parts.s1 + "\"position\":" + tokens[tokens.length - 1].substring(0, tokens[tokens.length - 1].indexOf("(")) + ",\"subject\":\""
					+ tokens[tokens.length - 1].substring(tokens[tokens.length - 1].indexOf("(") + 1, tokens[tokens.length - 1].indexOf(":")) + "\",\"query\":\""
					+ tokens[tokens.length - 1].substring(tokens[tokens.length - 1].indexOf(":") + 1, tokens[tokens.length - 1].indexOf(")")) + "\"";
			parts.s1 = parts.s1 + "}";

			parts.s2 = parts.s2 + "{\"begin\":" + tokens[tokens.length - 1].substring(0, tokens[tokens.length - 1].indexOf("(")) + ",\"end\":"
					+ (Integer.parseInt(tokens[tokens.length - 1].substring(0, tokens[tokens.length - 1].indexOf("("))) + offset);
			parts.s2 = parts.s2 + "}";
		}
		// }
		return parts;
	}

	public static String jsonq(String source) {
		String part = "";
		String tokens[];
		if (source.length() > 0) {
			tokens = source.split("[ ]+");
			if (tokens.length >= 1) {
				for (int z = 0; z < tokens.length - 1; z++) {
					part = part + "{";
					part = part + "\"position\":" + tokens[z].substring(0, tokens[z].indexOf("(")) + ",\"subject\":\"" + tokens[z].charAt(tokens[z].indexOf("(") + 1) + "\",\"query\":\""
					//		+ tokens[z].charAt(tokens[z].indexOf(":") + 1) + "\"";
							+ tokens[z].substring(tokens[z].indexOf(":") + 1) + "\"";
					part = part + "},";
				}
				part = part + "{";
				part = part + "\"position\":" + tokens[tokens.length - 1].substring(0, tokens[tokens.length - 1].indexOf("(")) + ",\"subject\":\""
						+ tokens[tokens.length - 1].charAt(tokens[tokens.length - 1].indexOf("(") + 1) + "\",\"query\":\""
				//		+ tokens[tokens.length - 1].charAt(tokens[tokens.length - 1].indexOf(":") + 1) + "\"";
						+ tokens[tokens.length - 1].substring(tokens[tokens.length - 1].indexOf(":") + 1) + "\"";
				part = part + "}";
				
				//substring für triplett
			}
		}
		return part;
	}

	// json-File generieren für Grafik.
	// {"mutations":[],"insertions":[],"deletions":[],"length":"707","start":"462","end":"556","found":true}
	// {"mutations":[{"position":66,"subject":"T","query":"-"},{"position":91,"subject":"-","query":"K"},{"position":92,"subject":"-","query":"G"},{"position":93,"subject":"-","query":"K"},{"position":178,"subject":"R","query":"L"},{"position":240,"subject":"E","query":"G"},{"position":265,"subject":"T","query":"A"},{"position":281,"subject":"Q","query":"R"}],"insertions":[{"begin":91,"end":93}],"deletions":[{"begin":66,"end":66}],"length":"376","start":"1","end":"376","found":true}
	public void jsongen(String filename) {

		MString all_ins = new MString();
		MString all_dels = new MString();
		MString all_sfs = new MString();
		MString all_dfs = new MString();
		MString all_stp = new MString();
		String muts = "";
		String m = "";
		String p = "";
		String ins = "";
		String dels = "";
		String sframe = "";
		String dframe = "";
		String stp = "";
		String wb = "";

		// json-parser-parts START
		// inserts
		all_ins = jsonp(mut.i_pos);
		// deletions
		all_dels = jsonp(mut.d_pos);
		// mutations
		m = jsonq(mut.r_pos);
		// positives
		p = jsonq(mut.p_pos);
		// fs-Zeug...
		all_sfs = jsonp(mut.sf_pos);
		all_dfs = jsonp(mut.df_pos);

		// stop-codons
		all_stp = jsonp(mut.stp_pos);

		// json-parser-parts ENDE

		String l = "\"length\":";
		l = l + "\"" + len1 + "\"";
		String s = "\"start\":\"" + startrow + "\"";
		String e = "\"end\":\"" + (bestrow - 2) + "\"";
		String f = "\"found\":";
		if (bestrow > startrow) {
			f = f + "true";
		} else {
			f = f + "false";
		}
		ins = "\"insertions\":[" + all_ins.s2 + "]";
		dels = "\"deletions\":[" + all_dels.s2 + "]";

		// all mutations zusammenbasteln
		muts = "\"mutations\":[" + m;
		if (p.length() > 0) {
			if (m.length() > 0) {
				muts = muts + "," + p;
			} else {
				muts = muts + p;
			}
		}
		/*
		 * Test für vereinfachtes JSON-File! if (all_ins.s1.length() > 0) { muts
		 * = muts + "," + all_ins.s1; } if (all_dels.s1.length() > 0) { muts =
		 * muts + "," + all_dels.s1; } if (all_sfs.s1.length() > 0) { muts =
		 * muts + "," + all_sfs.s1; } if (all_dfs.s1.length() > 0) { muts = muts
		 * + "," + all_dfs.s1; } if (all_stp.s1.length() > 0) { muts = muts +
		 * "," + all_stp.s1; }
		 */

		muts = muts + "]";

		sframe = "\"single frame-shifts\":[" + all_sfs.s2 + "]";
		dframe = "\"double frame-shifts\":[" + all_dfs.s2 + "]";
		stp = "\"stop-codons\":[" + all_stp.s2 + "]";

		// wilb-bases erfassen

		wb = "\"WB\":[";
		for (int i = 0; i < wb_counter.length; i++) {
			if (wb_counter[i] > 0) {
				switch (i) {
				case 0:
					wb = wb + "{\"K\":" + wb_counter[0] + "}";
					break;
				case 1:
					wb = wb + "{\"S\":" + wb_counter[1] + "}";
					break;
				case 2:
					wb = wb + "{\"Y\":" + wb_counter[2] + "}";
					break;
				case 3:
					wb = wb + "{\"M\":" + wb_counter[3] + "}";
					break;
				case 4:
					wb = wb + "{\"W\":" + wb_counter[4] + "}";
					break;
				case 5:
					wb = wb + "{\"R\":" + wb_counter[5] + "}";
					break;
				case 6:
					wb = wb + "{\"B\":" + wb_counter[6] + "}";
					break;
				case 7:
					wb = wb + "{\"D\":" + wb_counter[7] + "}";
					break;
				case 8:
					wb = wb + "{\"H\":" + wb_counter[8] + "}";
					break;
				case 9:
					wb = wb + "{\"V\":" + wb_counter[9] + "}";
					break;
				case 10:
					wb = wb + "{\"N\":" + wb_counter[10] + "}";
					break;
				}
				wb = wb + ",";
			}
		}
		if (wb.charAt(wb.length() - 1) == ',') {
			wb = wb.substring(0, wb.length() - 1);
		}
		wb = wb + "]";
		// daten in json-File schreiben

		Util.fill_file(filename, "{" + muts + "," + ins + "," + dels + "," + sframe + "," + dframe + "," + stp + "," + wb + "," + l + "," + s + "," + e + "," + f + "}");
	}

	// csv-gen für scalable grafik - noch dirty hack, aber tut schon mal ...

	public static void main(String[] args) throws IOException {

		int gop = -10;
		int gep = -1;
		int sfsp = -20;
		int dfsp = -30; // -40
		int stp = -4;
		boolean f1 = false;
		boolean noaffine = false;
		String seq1 = "";
		String seq2 = "";
		String n1 = "";
		String n2 = "";
		String json = "";
		String matrix = "BLOSUM62";
		String special = "";
		boolean worstwild = false;

		for (int i = 0; i < args.length; i++) {
			if (args[i].contains(".fasta")) {
				if (!f1) {
					n1 = readFile(args[i])[0];
					seq1 = readFile(args[i])[1];
					f1 = true;
				} else {
					n2 = readFile(args[i])[0];
					seq2 = readFile(args[i])[1];
				}
			}
			if (args[i].contains("g=")) {
				gop = Integer.parseInt(args[i].substring(2));
			}
			if (args[i].contains("e=")) {
				gep = Integer.parseInt(args[i].substring(2));
			}
			if (args[i].contains("s=")) {
				sfsp = Integer.parseInt(args[i].substring(2));
			}
			if (args[i].contains("d=")) {
				dfsp = Integer.parseInt(args[i].substring(2));
			}
			if (args[i].contains("p=")) {
				stp = Integer.parseInt(args[i].substring(2));
			}
			if (args[i].contains("m=")) {
				matrix = args[i].substring(2);
			}
			if (args[i].equals("na")) {
				noaffine = true;
				special = special + "|no affine gap costs|";
			}
			if (args[i].equals("w")) {
				worstwild = true;
				special = special + "|wild base mutations|";
			}
			if (args[i].equals("ng")) {
				gop = Integer.MIN_VALUE / 2;
				gep = Integer.MIN_VALUE / 2;
				special = special + "|no gaps at all|";
			}
			if (args[i].equals("nf")) {
				sfsp = Integer.MIN_VALUE / 2;
				dfsp = Integer.MIN_VALUE / 2;
				special = special + "|no frameshifts|";
			}
			if (args[i].contains("o=")) {
				json = args[i].substring(2);
			}
			if (noaffine) {
				gep = gop;
			}
		}

		SWAT sw = new SWAT(n1, seq1, n2, seq2, gop, gep, sfsp, dfsp, stp, matrix, special, worstwild);

		sw.align_wb();
		sw.setSize();
		sw.setMut();
		sw.setStop();
		System.out.println("SWAT-Alignment-Protocol:");
		System.out.println("");
		sw.showNames();
		System.out.println("");
		sw.showSLength();
		System.out.println("");
		System.out.println("");
		sw.showParameters();
		System.out.println("");
		System.out.println("");
		sw.showScore();
		sw.showSize();
		sw.showALength();
		System.out.println("");
		System.out.println("");
		sw.showAlignment();
		System.out.println("");
		System.out.println("");
		sw.showIdens();
		sw.showPosis();
		sw.showGaps();
		System.out.println("");
		System.out.println("");
		sw.showMut();
		System.out.println("");
		sw.showStop();
		System.out.println("");
		sw.showWBs();
		if (json.length() > 0) {
			sw.jsongen(json);
		}
		System.out.println("");

	}
}
