import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

public class Util {
	public static int max(int a, int b, int c, int d) {
		int val = d;
		if (a > val)
			val = a;
		if (b > val)
			val = b;
		if (c > val)
			val = c;

		return val;
	}

	public static void fill_file(String filename, String Dataset) {
		try {
			BufferedWriter s_data = new BufferedWriter(new FileWriter(filename, true));

			s_data.flush();
			s_data.write(Dataset);
			s_data.newLine();
			// System.out.println(line);

			s_data.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public static String toCsv(String[][] arr) {
		String csvString = "";
		String allString = "";
		for (int j = 0; j < arr.length; j++) {
			csvString = "" + arr[j][0];
			for (int i = 1; i < arr[j].length; i++) {
				csvString = csvString + "," + arr[j][i];
			}
			csvString += "\n";
			allString = allString + csvString;
		}
		return allString;
	}

}
