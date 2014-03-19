import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

public class Matrix {

	static int[][] mgen(String mname, int stp) throws IOException {
		BufferedReader in = new BufferedReader(new FileReader(mname));
		int[][] matrix = new int[26][26];
		int j = 0, k = 0;
		String zeile = "";
		zeile = in.readLine();
		while ((zeile = in.readLine()) != null) {

			String[] tokens = zeile.split("[ ]+");

			for (int i = 0; i < tokens.length; i++) {
				try {
					matrix[k][j] = Integer.parseInt(tokens[i]);
					j++;
				} catch (Exception e) { // Pokemon....
				}
			}
			j=0;
			k++;
		}
		in.close();
		for(int i=0;i<=25;i++){
			matrix[i][24]=-100000;
			matrix[i][25]= stp;
			matrix[24][i]=-100000;
			matrix[25][i]= stp;
		}
		matrix[24][24]=-100000;
		matrix[25][25]= 1; // was wenn Stp auf Stp trifft?
	
		return matrix;
	}
}
