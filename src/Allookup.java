public class Allookup {

	public static int lookup(char letter) {
		//  A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V, B, Z, X, *, |, .
		// 00,01,02,03,04,05,06,07,08,09,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25
		switch (letter) {
		case 'A':
			return 0;
		case 'R':
			return 1;
		case 'N':
			return 2;
		case 'D':
			return 3;
		case 'C':
			return 4;
		case 'Q':
			return 5;
		case 'E':
			return 6;
		case 'G':
			return 7;
		case 'H':
			return 8;
		case 'I':
			return 9;
		case 'L':
			return 10;
		case 'K':
			return 11;
		case 'M':
			return 12;
		case 'F':
			return 13;
		case 'P':
			return 14;
		case 'S':
			return 15;
		case 'T':
			return 16;
		case 'W':
			return 17;
		case 'Y':
			return 18;
		case 'V':
			return 19;
		case 'B':
			return 20;
		case 'Z':
			return 21;
		case 'X':
			return 22;
		case '*':
			return 23;
		case '|':
			return 24;
		case '.':
			return 25;

		default:
			return -1;
		}

	}

	public static int lookup2(char letter) {

		switch (letter) {
		case 'A':
			return 0;
		case 'T':
			return 1;
		case 'G':
			return 2;
		case 'C':
			return 3;
		case '|':
			return 4;
		default:
			return -1;
		}
	}
	
	public static String wb_lookup(char wb){
		switch(wb){
		case 'A':
		case 'C':
		case 'G':
		case 'T':
		case 'U': return ""+wb;
		case 'K': return "GU";
		case 'S': return "GC";
		case 'Y': return "UC"; 
		case 'M': return "AC";
		case 'W': return "AU";
		case 'R': return "GA";
		case 'B': return "GUC";
		case 'D': return "GAU";
		case 'H': return "ACU";
		case 'V': return "GCA";
		case 'N': return "AGCU";
		default: return "";
		}
	}

	public static char trilookup(String triplet) {
		switch (triplet) {
		case "UUC":
		case "TTC":
		case "UUU":
		case "TTT":
			return 'F';

		case "UUA":
		case "UUG":
		case "CUU":
		case "CUA":
		case "CUC":
		case "CUG":
		case "TTA":
		case "TTG":
		case "CTT":
		case "CTA":
		case "CTC":
		case "CTG":
			return 'L';

		case "UCU":
		case "UCA":
		case "UCG":
		case "UCC":
		case "AGU":
		case "AGC":
		case "TCT":
		case "TCA":
		case "TCG":
		case "TCC":
		case "AGT":
			return 'S';

		case "UAU":
		case "UAC":
		case "TAT":
		case "TAC":
			return 'Y';

		case "UGU":
		case "UGC":
		case "TGT":
		case "TGC":
			return 'C';

		case "UAA":
		case "UGA":
		case "UAG":
		case "TAA":
		case "TGA":
		case "TAG":
			return '*'; // FIX FÜR STOPCODON!!!!

		case "UGG":
		case "TGG":
			return 'W';

		case "CCU":
		case "CCC":
		case "CCG":
		case "CCA":
		case "CCT":
			return 'P';

		case "CAU":
		case "CAC":
		case "CAT":
			return 'H';

		case "CGU":
		case "CGC":
		case "CGA":
		case "CGG":
		case "AGG":
		case "AGA":
		case "CGT":
			return 'R';

		case "CAA":
		case "CAG":
			return 'Q';

		case "AUU":
		case "AUC":
		case "AUA":
		case "ATT":
		case "ATC":
		case "ATA":
			return 'I';

		case "ACU":
		case "ACC":
		case "ACA":
		case "ACG":
		case "ACT":
			return 'T';

		case "AAU":
		case "AAC":
		case "AAT":
			return 'N';

		case "AAA":
		case "AAG":
			return 'K';

		case "AUG":
		case "ATG":
			return 'M';

		case "GUU":
		case "GUC":
		case "GUG":
		case "GUA":
		case "GTT":
		case "GTC":
		case "GTG":
		case "GTA":
			return 'V';

		case "GCU":
		case "GCC":
		case "GCA":
		case "GCG":
		case "GCT":
			return 'A';

		case "GAU":
		case "GAC":
		case "GAT":
			return 'D';

		case "GGU":
		case "GGC":
		case "GGA":
		case "GGG":
		case "GGT":
			return 'G';

		case "GAA":
		case "GAG":
			return 'E';
		case "---":
			return '_';
		default: {
//			System.out.println("fail!");
			return 'X'; // mal 'X' statt '|' für submitted wb.....

		}
		}
	}

	public static String p2pro(char c){
		switch(c){
		case 'F' : return "Phe";
		case 'I' : return "Ile";
		case 'P' : return "Pro";
		case 'T' : return "Thr";
		case 'A' : return "Ala";
		case 'V' : return "Val";
		case 'D' : return "Asp";
		case 'E' : return "Glu";
		case 'G' : return "Gly";
		case 'R' : return "Arg";
		case 'S' : return "Ser";
		case 'N' : return "Asn";
		case 'Q' : return "Gln";
		case 'H' : return "His";
		case 'Y' : return "Tyr";
		case 'C' : return "Cys";
		case 'W' : return "Trp";
		case 'M' : return "Met";
		case 'L' : return "Leu";
		case 'K' : return "Lys";
		case '.' : return "STP";
		case '_' : return "---";
		default  : return "xxx";
		}
	}
}
