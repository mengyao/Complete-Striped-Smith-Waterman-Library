package ssw;

public class Example {
	public static void main(String[] args) {
		System.loadLibrary("sswjni");
		/*
		System.out.println("Creating profile");
		long ptr = Aligner.initprofile(new byte[] {0,1,2,0}, new byte[] {
				1,0,0,
				0,1,0,
				0,0,1}, 3, 0);
		System.out.println("Aligning with profile " + Long.toHexString(ptr));
		Alignment aln = Aligner.align(ptr, new byte[] { 0, 0, 1, 2, 2, 2, 2 }, (byte)10, (byte)1, (byte)0, (short)0, 0, 15);
		
		System.out.println("Alignment " + aln.toString());
		System.out.println(String.format("score1=%d score2=%d, ref=%d,%d read=%d,%d refend2=%d",
				aln.score1, aln.score2, aln.ref_begin1, aln.ref_end1, aln.read_begin1, aln.read_end1, aln.ref_end2));
		Aligner.destroyprofile(ptr);
		*/
		int[][] score = new int[128][128];
		for (int i = 0; i < 128; i++) {
			for (int j = 0; j < 128; j++) {
				if (i == j) score[i][j] = 2;
				else score[i][j] = -2;
			}
		}
		System.out.println("Aligning nucleotides");
		Alignment aln = Aligner.align("CTGAGCCGGTAAATC".getBytes(), "CAGCCTTTCTGACCCGGAAATCAAAATAGGCACAACAAA".getBytes(), score, 3, 1, true);
		if (aln == null) {
			throw new RuntimeException();
		}
		System.out.print(String.format("score1=%d ", aln.score1));
		System.out.print(String.format("score2=%d ", aln.score2));
		System.out.print(String.format("ref_begin1=%d ", aln.ref_begin1));
		System.out.print(String.format("ref_end1=%d ", aln.ref_end1));
		System.out.print(String.format("read_begin1=%d ", aln.read_begin1));
		System.out.print(String.format("read_end1=%d ", aln.read_end1));
		System.out.print(String.format("ref_end2=%d ", aln.ref_end2));
		if (aln.cigar != null) System.out.print(aln.cigar);
	}
}
