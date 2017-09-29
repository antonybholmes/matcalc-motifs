package edu.columbia.rdf.matcalc.toolbox.motifs;

import java.io.IOException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.jebtk.bioinformatics.dna.Sequence;
import org.jebtk.bioinformatics.genomic.Gene;
import org.jebtk.bioinformatics.genomic.GenesDb;
import org.jebtk.bioinformatics.genomic.GenomicRegion;
import org.jebtk.bioinformatics.genomic.SequenceRegion;
import org.jebtk.bioinformatics.motifs.Motif;
import org.jebtk.core.Mathematics;
import org.jebtk.core.collections.CollectionUtils;
import org.jebtk.core.collections.DefaultHashMap;
import org.jebtk.core.collections.DefaultHashMapCreator;
import org.jebtk.core.collections.TreeSetCreator;
import org.jebtk.core.io.PathUtils;
import org.jebtk.bioinformatics.ui.groups.Group;




public class MotifSearch {
	private static final double ROUNDING_ERROR = 1E-10;

	public static final int BUFFER_SIZE = 100000;

	public static final double[] FG = new double[BUFFER_SIZE];

	public static final double[] LLKRF = new double[BUFFER_SIZE];
	public static final double[] LLKRR = new double[BUFFER_SIZE];

	public static void main(String[] args) throws IOException {
		//albertoTest();
		check1Test();
	}

	public static void albertoTest() throws IOException {

		/*
		List<Sequence> foregroundSequences = 
				SequenceRegion.parseFasta(PathUtils.getPath("/ifs/scratch/cancer/Lab_RDF/abh2138/matlab/top_250.fasta"));

		List<Sequence> foregroundRevCompSeqs =
				Sequence.reverseComplement(foregroundSequences);

		List<Sequence> backgroundSequences = 
				Sequence.parseFasta(PathUtils.getPath("/ifs/scratch/cancer/Lab_RDF/abh2138/matlab/bottom_250.fasta"));

		List<Sequence> backgroundRevCompSeqs = 
				Sequence.reverseComplement(backgroundSequences);

		List<Motif> motifs = new ArrayList<Motif>(); 

		//motifs.addAll(Motif.parseMotifs(Resources.getGZIPReader("res/motifs.txt.gz")));
		motifs.addAll(Motif.parseMotifs(PathUtils.getPath("/ifs/scratch/cancer/Lab_RDF/abh2138/matlab/stat1.motif"), "test"));
		motifs.addAll(Motif.parseMotifs(PathUtils.getPath("/ifs/scratch/cancer/Lab_RDF/abh2138/matlab/ebf1.motif"), "test"));
		motifs.addAll(Motif.parseMotifs(PathUtils.getPath("/ifs/scratch/cancer/Lab_RDF/abh2138/matlab/bcl6.motif"), "test"));


		// Find the best score for each motif in each sequence


		boolean[] goldStandard = 
				createGoldStandard(foregroundSequences.size(), backgroundSequences.size());

		int threshold = 0;

		double minSensitivity = 0.3;
		double minSpecificity = 0.5;

		int ns = foregroundSequences.size() + 
				backgroundSequences.size();

		double[] bestScores = new double[ns];

		int c = 0;


		for (Motif m : motifs) {
			System.err.println(m.getName());

			int w = m.getBaseCount();

			bestScores(m,
					w,
					foregroundSequences,
					foregroundRevCompSeqs,
					threshold,
					bestScores, 
					0);


			bestScores(m,
					w, 
					backgroundSequences, 
					backgroundRevCompSeqs,
					threshold,
					bestScores,
					foregroundSequences.size());

			Stats stats = enrichmentByMinError(bestScores, 
					goldStandard, 
					minSensitivity, 
					minSpecificity);

			LOG.debug("Min error: {}", stats.error);

			// Once the error has been minimized, we can calculate the rest of the
			// stats
			if (stats.error < 1) {
				double ppv = (double)stats.truePositive / (double)(stats.truePositive + stats.falsePositive);

				Hypergeometric hyg = new Hypergeometric();

				LOG.debug("tp: {}, tn: {}, fp: {}, fn: {}, n: {}", 
						stats.truePositive,
						stats.trueNegative, 
						stats.falsePositive,
						stats.falseNegative, 
						ns);

				//double p = hyg.cdfTwoTail(stats.truePositive, ns / 2, stats.truePositive + stats.falsePositive, ns);
				double p = hyg.cdfOneTail(stats.truePositive, ns / 2, stats.truePositive + stats.falsePositive, ns);

				LOG.debug("HygP: {}", p);

				if (p <= 0.01) {
					++c;
				}
			}
		}

		System.err.println("all the p " + c);
		 */


		/*

		// RBPJ motif
		Motif motif = motifs.get(1347); //motifs.get(541); //Motif.parseJaspar(PathUtils.getPath("rbpj.jaspar")).get(0);

		threshold = 12;

		// See how prevalent the motif is in all the promoters

		double motifBackground = 0;

		int w = motif.getBaseCount();

		double bgpwm = getBgPWMScore(w);

		Sequence reverseComplement;

		for (Promoter p : promoters) {
			int n = p.getLength();

			List<Double> bgscores = Mathematics.repeat(bgpwm, n);

			reverseComplement = Sequence.reverseComplement(p);

			List<Double> scores = search(p, reverseComplement, motif, bgscores, n, w);

			//Sequence scoredSequence = showScores(p, scores, threshold);

			List<PotentialBindingSite> sites = 
					getPotentialBindingSites(p, scores);

			List<PotentialBindingSite> thresholdedSites = 
					thresholdSites(sites, threshold);

			if (thresholdedSites.size() > 0) {
				//System.err.println(p.getName());

				++motifBackground;
			}
		}

		double p = motifBackground / promoters.size();

		System.err.println("Background p " + motifBackground + "/" + promoters.size() + " " + p);

		System.exit(0);


		//motifs = ArrayUtils.toList(motif);

		System.err.println(motif.getId());

		Sequence sequence = promoters.get(13962);

		reverseComplement = Sequence.reverseComplement(sequence);

		System.err.println("hasdasd " + sequence.getName());

		int n = sequence.getLength();

		int c = 0;

		List<Double> bgscores = Mathematics.repeat(bgpwm, n);

		for (Motif m : motifs) {
			w = motif.getBaseCount();

			bgpwm = Math.pow(0.25, w);

			List<Double> scores = search(sequence, reverseComplement, m, bgscores, n, w);

			Sequence scoredSequence = showScores(sequence, scores, threshold);

			List<PotentialBindingSite> sites = 
					getPotentialBindingSites(sequence, scores);

			List<PotentialBindingSite> thresholdedSites = 
					thresholdSites(sites, threshold);

			if (thresholdedSites.size() > 0) {
				System.err.println(m.getName());

				//PotentialBindingSite best = getBestBindingSites(sites);



				//System.err.println(scores.toString());
				System.err.println(scoredSequence.getBases());
				//System.err.println("n:" + sites.size());
				System.err.println(thresholdedSites.toString());
				//System.err.println(best);

				++c;
			}
		}

		System.err.println("c " + c + " / " + motifs.size());
		 */
	}

	public static void check1Test() throws IOException {

		//List<Sequence> foregroundSequences =
		//		Sequence.parseFasta(PathUtils.getPath("test2.fasta"));

		//List<Sequence> foregroundSequences =
		//		Sequence.parseFasta(PathUtils.getPath("all_sorted_foreground.fasta"));

		List<Sequence> foregroundSequences = 
				Sequence.parseFasta(PathUtils.getPath("chek1.fasta"));

		//List<Sequence> foregroundSequences = 
		//		Sequence.parseFasta(PathUtils.getPath("foreground.fasta"));

		//List<Sequence> foregroundRevCompSeqs =
		//		Sequence.reverseComplement(foregroundSequences);

		//SequenceRegion.writeFasta(backgroundSequences, PathUtils.getPath("background.fasta"));

		//List<Sequence> foregroundSequences =
		//		Sequence.parseFasta(PathUtils.getPath("foreground.fasta"));

		//List<Sequence> backgroundSequences =
		//		Sequence.parseFasta(PathUtils.getPath("background.fasta"));

		//System.exit(0);

		//List<Motif> motifs = 
		//		Motif.parseMotifs(PathUtils.getPath("/ifs/scratch/cancer/Lab_RDF/abh2138/matlab/met4.motif"));

		List<Motif> motifs = new ArrayList<Motif>(); 

		//motifs.addAll(Motif.parseMotifs(Resources.getGZIPReader("res/motifs.txt.gz")));
		//motifs.addAll(Motif.parseMotifs(PathUtils.getPath("res/motifs/bcl6_m0.motif")));
		//motifs.addAll(Motif.parseMotifs(PathUtils.getPath("res/motifs/bcl6_m00424.motif")));
		//motifs.addAll(Motif.parseMotifs(PathUtils.getPath("res/motifs/bcl6_m2.motif")));
		motifs.addAll(Motif.parseMotifs(PathUtils.getPath("res/motifs/bcl6_b6bs.motif"), "test"));

		int n = foregroundSequences.get(0).getLength();

		//double[] llkrf = new double[n];
		//double[] llkrr = new double[n];

		double threshold = 0;

		for (Motif motif : motifs) {
			//System.err.println(m.getName());

			int w = motif.getBaseCount();

			double bgscore = motif.getBgPwm();
			//double[] bgscores = Mathematics.repeatArray(m.getBgPwm(), n);
			
			double[][] pwm = motif.getPwm();

			for (int i = 0; i < foregroundSequences.size(); ++i) {
				char[] seq = foregroundSequences.get(i).getChars();
				char[] revCompSeq = foregroundSequences.get(i).getChars();

				byte[] iSeq = Sequence.seqToIndexSeq(seq);
				byte[] iRevCompSeq = Sequence.seqToIndexSeq(revCompSeq);

				//System.err.println(sequence);

				MotifSearch.search(iSeq, 
						iRevCompSeq,
						seq.length,
						pwm, 
						bgscore,
						w, 
						threshold,
						LLKRF, 
						LLKRR);


				System.err.println("nz: " + Mathematics.nonZeroIndices(LLKRR));

				//Sequence showSequence = MotifSearch.showScores(sequence, 
				///		llkrf, 
				//		llkrr,
				//		0);

				List<BindingSite> sites = 
						MotifSearch.getBindingSites(seq, n, w, LLKRF, LLKRR);

				for (BindingSite site : sites) {
					Sequence showSequence = MotifSearch.showScores(seq, 
							site,
							n);

					System.err.println("offset " + site.getOffset());
				}
			}
		}
	}

	/**
	 * To generate max scores we must look for the highest value in a 
	 * block of zero values immediately transitioning to a run of scores
	 * equal to the length of the motif. Scores must have been thresholded
	 * and log transformed to ensure that small values are floored to zero.
	 * 
	 * @param scores		The scores.
	 * @param n				The sequence length.
	 * @param w				The width of the motif.
	 * @return				The max score found for the motif, or 0 if
	 * 						nothing was found.
	 */
	private static double max(final double[] scores, int n, int w) {
		double max = 0;

		int i = 0;
		int wi = w - 1;

		int c = n - w + 1;

		while (i < c) {
			/*
			if (previous <= 0.0 && scores[i] > 0) {
				// test the length

				int s = (int)scores[i];

				boolean found = true;

				for (int j = i; j < i + w; ++j) {
					if ((int)scores[j] < s) {
						found = false;
						break;
					}
				}

				if (found) {
					if (scores[i] > max) {
						max = scores[i];
						System.err.println("max " + i + " " + max);

					}

					// There must be at least 1 one base gap between
					// one binding site and the next.
					i += wi;
				}
			}
			 */

			if (scores[i] > 0) {
				// test the length

				if (scores[i] > max) {
					max = scores[i];
				}

				// There must be at least 1 one base gap between
				// one binding site and the next.
				i += wi;
			}

			//previous = scores[i];

			++i;
		}

		return max;
	}

	public static Stats enrichmentByMinError(final double[] bestScores, 
			final boolean[] goldStandard,
			double minSensitivity,
			double minSpecificity) {

		List<Double> thresholds = 
				CollectionUtils.reverse(CollectionUtils.sort(CollectionUtils.unique(Mathematics.gtZero(Mathematics.floor(CollectionUtils.head(bestScores, bestScores.length / 2))))));

		//LOG.debug("Min Sensitivity: {}, Min Specificity: {}", 
		//		minSensitivity, 
		//		minSpecificity);

		//LOG.debug("Thresholds: {}", thresholds);

		int numSequences = bestScores.length;

		Stats stats = new Stats();

		for (double threshold : thresholds) {
			//LOG.debug("Threshold: {}", threshold);

			int tp = 0;
			int fp = 0;
			int fn = 0;
			int tn = 0;

			for (int i = 0; i < numSequences; ++i) {
				boolean aboveThreshold = bestScores[i] >= threshold;
				boolean isForeground = goldStandard[i]; //goldStandard[i] == 1;

				if (aboveThreshold && isForeground) {
					++tp;
				}

				if (aboveThreshold && !isForeground) {
					++fp;
				}

				if (!aboveThreshold && isForeground) {
					++fn;
				}

				if (!aboveThreshold && !isForeground) {
					++tn;
				}
			}

			double sensitivity = tp / (double)(tp + fn);
			double specificity = tn / (double)(tn + fp);

			double error = 1.0 - (sensitivity + specificity) / 2.0;
			//double error = 1.0 - sensitivity;
			//double error = 1.0 - specificity;

			//MotifsModule.LOG.info("Sensitivity: {}, Specificity: {}, Error: {}, Threshold: {}", 
			//		sensitivity, 
			//		specificity, 
			//		error,
			//		threshold);
								
			//LOG.debug("tp: {}, fp: {}, fn: {}, tn: {}", tp, fp, fn, tn);

			if (sensitivity >= minSensitivity && specificity >= minSpecificity) {
				if (error <= stats.error) {
					stats.sensitivity = sensitivity;
					stats.specificity = specificity;
					stats.error = error;
					stats.truePositive = tp;
					stats.trueNegative = tn;
					stats.falsePositive = fp;
					stats.falseNegative = fn;
				}
			}
		}

		return stats;
	}

	public static List<Stats> enrichment(final double[] bestScores, 
			final boolean[] goldStandard) {

		List<Double> thresholds = 
				CollectionUtils.reverse(CollectionUtils.sort(CollectionUtils.unique(Mathematics.gtZero(Mathematics.floor(CollectionUtils.head(bestScores, bestScores.length / 2))))));

		//LOG.debug("Min Sensitivity: {}, Min Specificity: {}", 
		//		minSensitivity, 
		//		minSpecificity);

		//LOG.debug("Thresholds: {}", thresholds);

		int numSequences = bestScores.length;

		List<Stats> allStats = new ArrayList<Stats>();

		for (double threshold : thresholds) {
			Stats stats = new Stats();

			allStats.add(stats);

			int tp = 0;
			int fp = 0;
			int fn = 0;
			int tn = 0;

			for (int i = 0; i < numSequences; ++i) {
				boolean aboveThreshold = bestScores[i] >= threshold;
				boolean isForeground = goldStandard[i]; //goldStandard[i] == 1;

				if (aboveThreshold && isForeground) {
					++tp;
				}

				if (aboveThreshold && !isForeground) {
					++fp;
				}

				if (!aboveThreshold && isForeground) {
					++fn;
				}

				if (!aboveThreshold && !isForeground) {
					++tn;
				}
			}

			double sensitivity = tp / (double)(tp + fn);
			double specificity = tn / (double)(tn + fp);

			double error = 1.0 - (sensitivity + specificity) / 2.0;

			/*
			LOG.debug("Sensitivity: {}, Specificity: {}, Error: {}, Threshold: {}", 
					sensitivity, 
					specificity, 
					error,
					threshold);

			LOG.debug("tp: {}, fp: {}, fn: {}, tn: {}", tp, fp, fn, tn);
			 */

			stats.threshold = threshold;
			stats.sensitivity = sensitivity;
			stats.specificity = specificity;
			stats.error = error;
			stats.truePositive = tp;
			stats.trueNegative = tn;
			stats.falsePositive = fp;
			stats.falseNegative = fn;
		}

		return allStats;
	}

	/**
	 * The gold standard indicates when we should see true positives in
	 * the concatenated top and bottom sequences. It is simply 1 for the
	 * top n sequences and 0 for the bottom n sequences.
	 * 
	 * @param foregroundSize
	 * @param backgroundSize
	 * @return
	 */
	public static boolean[] createGoldStandard(int foregroundSize, 
			int backgroundSize) {
		boolean[] ret = new boolean[foregroundSize + backgroundSize];

		CollectionUtils.copyValue(true, ret, foregroundSize, 0);
		CollectionUtils.copyValue(false, ret, backgroundSize, foregroundSize);

		return ret;
	}

	public static <X extends SequenceRegion> void bestScores(final Motif motif,
			int w, 
			final List<X> seqs,
			final List<X> revCompSeqs,
			double threshold,
			double[] bestScores) {
		bestScores(motif, w, seqs, revCompSeqs, threshold, 0, bestScores);
	}

	/**
	 * 
	 * @param <X>
	 * @param motif
	 * @param n
	 * @param w
	 * @param seqs
	 * @param revCompSeqs
	 * @param bestScores
	 * @param offset				Offset to start writing in best scores.
	 */
	public static <X extends SequenceRegion> void bestScores(final Motif motif,
			int w, 
			final List<X> seqs,
			final List<X> revCompSeqs,
			double threshold,
			int offset,
			double[] bestScores) {

		double bgscore = motif.getBgPwm();
		
		double[][] pwm = motif.getPwm();

		//double[] llkrf = new double[BUFFER_SIZE];
		//double[] llkrr = new double[BUFFER_SIZE];


		byte[][] iSeqs = SequenceRegion.toIndex(seqs);
		byte[][] iRevCompSeqs = SequenceRegion.toIndex(revCompSeqs);

		for (int i = 0; i < iSeqs.length; ++i) {
			byte[] iSeq = iSeqs[i];
			byte[] iRevCompSeq = iRevCompSeqs[i];

			int n = iSeq.length;

			search(iSeq, iRevCompSeq, n, pwm, bgscore, w, threshold, LLKRF, LLKRR);

			double max = Math.max(max(LLKRF, n, w), max(LLKRR, n, w));

			bestScores[i + offset] = max;
		}
	}
	
	public static <X extends SequenceRegion> void bestScores(final Motif motif,
			int w, 
			byte[][] iSeqs,
			byte[][] iRevCompSeqs,
			double threshold,
			int offset,
			double[] bestScores) {

		double bgscore = motif.getBgPwm();
		
		double[][] pwm = motif.getPwm();

		for (int i = 0; i < iSeqs.length; ++i) {
			byte[] iSeq = iSeqs[i];
			byte[] iRevCompSeq = iRevCompSeqs[i];

			int n = iSeq.length;

			search(iSeq, iRevCompSeq, n, pwm, bgscore, w, threshold, LLKRF, LLKRR);

			double max = Math.max(max(LLKRF, n, w), max(LLKRR, n, w));

			bestScores[i + offset] = max;
		}
	}
	
	public static <X extends SequenceRegion> void bestScores(final Motif motif,
			int w, 
			byte[][] iSeqs,
			Map<Integer, Map<Integer, Collection<Integer>>> tripletMap,
			byte[][] iRevCompSeqs,
			Map<Integer, Map<Integer, Collection<Integer>>> revTripletMap,
			double threshold,
			int offset,
			double[] bestScores) {

		double bgscore = motif.getBgPwm();
		
		double[][] pwm = motif.getPwm();

		for (int triplet : motif.getTriplets()) {
			if (tripletMap.containsKey(triplet)) {
				Map<Integer, Collection<Integer>> seqIndexMap = 
						tripletMap.get(triplet);

				for (int i : seqIndexMap.keySet()) {
					byte[] iSeq = iSeqs[i];
					
					Collection<Integer> startLocations = seqIndexMap.get(i);
					
					//System.err.println("trip " + triplet + " " + startLocations);

					int n = iSeq.length;

					search(iSeq, startLocations, n, pwm, bgscore, w, threshold, LLKRF);

					double max = max(LLKRF, n, w);

					bestScores[i + offset] = max;
				}
			}
		}

		for (int triplet : motif.getTriplets()) {
			if (revTripletMap.containsKey(triplet)) {
				Map<Integer, Collection<Integer>> seqIndexMap = 
						revTripletMap.get(triplet);

				for (int i : seqIndexMap.keySet()) {
					byte[] iRevCompSeq = iRevCompSeqs[i];
					
					Collection<Integer> startLocations = seqIndexMap.get(i);

					int n = iRevCompSeq.length;

					search(iRevCompSeq, startLocations, n, pwm, bgscore, w, threshold, LLKRR);

					// Need to reverse llkrr so it is in the forward direction, since
					// the buffer may exceed the sequence length, only reverse the
					// first n bases.
					CollectionUtils.reverse(LLKRR, 0, n);

					double max = max(LLKRR, n, w);

					bestScores[i + offset] = Math.max(bestScores[i + offset], max);
				}
			}
		}

		/*
		 * Older linear method
		 * 
		for (int i = 0; i < iSeqs.length; ++i) {
			byte[] iSeq = iSeqs[i];
			byte[] iRevCompSeq = iRevCompSeqs[i];

			int n = iSeq.length;

			//double[] llkrf = new double[n];
			//double[] llkrr = new double[n];

			search(iSeq, iRevCompSeq, n, motif, bgscore, w, threshold, LLKRF, LLKRR);

			double max = Math.max(max(LLKRF, n, w), max(LLKRR, n, w));

			bestScores[i + offset] = max;
		}
		 */
	}

	public static Map<Integer, Map<Integer, Set<Integer>>> mapTriplets(byte[][] seqs) {
		Map<Integer, Map<Integer, Set<Integer>>> ret =
				DefaultHashMap.create(new DefaultHashMapCreator<Integer, Set<Integer>>(new TreeSetCreator<Integer>()));

		for (int i = 0; i < seqs.length; ++i) {
			byte[] seq = seqs[i];

			int triplet = 0;

			for (int j = 0; j < seq.length - 3; ++j) {
				triplet = seq[j] * 100 + seq[j + 1] * 10 + seq[j + 2];

				// Index the triplet location in sequence i
				ret.get(triplet).get(i).add(j);
			}
		}

		return ret;
	}

	/**
	 * Get the best scores on one set of strands, e.g. if we only want
	 * to search the forward strands.
	 * 
	 * @param motif
	 * @param n
	 * @param w
	 * @param sequences
	 * @param bestScores
	 * @param offset
	 */
	public static void bestScores(final Motif motif,
			int w, 
			final List<Sequence> sequences,
			double threshold,
			int offset,
			double[] bestScores) {


		double bgscore = motif.getBgPwm(); //Mathematics.repeatArray(m.getBgPwm(), n);

		double[][] pwm = motif.getPwm();
		
		// Should hold most sequences we are searching for
		//double[] llkrf = new double[BUFFER_SIZE];

		for (int i = 0; i < sequences.size(); ++i) {
			char[] seq = sequences.get(i).toArray();
			byte[] iSeq = Sequence.seqToIndexSeq(seq);


			int n = seq.length;
			//double[] llkrf = new double[n];

			//System.err.println(sequence.getBases());
			//System.err.println(reverse.getBases());

			search(iSeq, n, pwm, bgscore, w, threshold, LLKRF);

			double max = Mathematics.max(LLKRF, 0, n);

			bestScores[i + offset] = max;
		}

		//System.err.println(bestScores);

		//return bestScores;
	}

	public static double bestScores(Motif motif,
			int w, 
			final byte[] sequence,
			final byte[] reverseComplement,
			int n,
			double threshold) {

		double bgscore = motif.getBgPwm(); //Mathematics.repeatArray(m.getBgPwm(), n);
		//double[] llkrf = new double[n];
		//double[] llkrr = new double[n];

		search(sequence,
				reverseComplement,
				n,
				motif.getPwm(), 
				bgscore,
				w,
				threshold,
				LLKRF, 
				LLKRR);

		double max = Math.max(max(LLKRF, n, w), max(LLKRR, n, w));

		return max;
	}

	/**
	 * Best scores along one strand.
	 * 
	 * @param motif
	 * @param n
	 * @param w
	 * @param sequence
	 * @param threshold
	 * @return
	 */
	public static double bestScores(Motif motif,
			int w, 
			final byte[] sequence,
			int n,
			double threshold) {

		double bgscore = motif.getBgPwm(); //Mathematics.repeatArray(m.getBgPwm(), n);
		//double[] llkrf = new double[n];

		search(sequence, n, motif.getPwm(), bgscore, w, threshold, LLKRF);

		double max = max(LLKRF, n, w);

		return max;
	}


	/**
	 * Best scores are determined from the int of the score of interest
	 * Thus the best score must be at least 1 to be effective.
	 * 
	 * @param bestScore
	 * @param llkrf
	 * @param llkrr
	 * @param w				the width of the motif.
	 * @return
	 */
	public static List<Integer> getBestScoreStartIndices(final char[] sequence,
			int n,
			double bestScore, 
			final double[] llkrf,
			final double[] llkrr,
			int w) {
		List<Integer> indices = new ArrayList<Integer>();

		double score = 0;

		int i = 0;

		while (i < n) {
			score = Math.max(llkrf[i], llkrr[i]);

			// Since we are dealing with double arithmetic, we
			// must test that the score is very close to the best
			// score, not that is equal to account for minor
			// rounding errors
			if (Math.abs(score - bestScore) <= ROUNDING_ERROR) {
				indices.add(i);

				// skip past this motif, there is no point
				// testing the other bases it occupies
				i += w;
			}

			// motifs must at least 1 bp apart hence w + 1
			++i;
		}

		return indices;
	}

	public static double getBgPWMScore(int w) {
		return Math.pow(0.25, w);
	}

	/**
	 * Returns the number of binding sites from the score array.
	 * 
	 * @param scores			The score for each base.
	 * @param w					The length of the motif.
	 * @param threshold			The threshold to consider a score significant.
	 * @return
	 */
	public static int getNumSites(List<Double> scores, 
			int w, 
			double threshold) {
		int c = 0;

		for (double score : scores) {
			if (Math.abs(score) > threshold) {
				++c;
			}
		}

		return c / w;	
	}

	/**
	 * Returns the list of binding sites in the given sequence from the scores
	 * at a particular threshold. Offsets are relative to the start of the
	 * sequence.
	 * 
	 * @param sequence		The sequence to search.
	 * @param n				Length of the sequence.
	 * @param llkrf			Scores on the forward strand.
	 * @param llkrr			Scores on the reverse strand.
	 * 
	 * @return
	 */
	public static List<BindingSite> getBindingSites(final char[] sequence,
			int n,
			int w,
			final double[] llkrf, 
			final double[] llkrr) {

		List<BindingSite> sites = processSite(sequence, n, w, '+', llkrf);

		// Append the reverse sequence matches
		sites.addAll(processSite(sequence, n, w, '-', llkrr));

		return sites;
	}

	/**
	 * Add binding sites to list.
	 * 
	 * @param sequence		The sequence being searched.
	 * @param w				The motif length
	 * @param strand
	 * @param scores
	 * @param sites
	 * @return 
	 */
	public static List<BindingSite> processSite(final char[] sequence,
			int n,
			int w,
			char strand,
			final double[] scores) {

		List<BindingSite> sites = new ArrayList<BindingSite>();

		for (int i = 0; i < n - w; ++i) {
			if (scores[i] > 0) {
				//BindingSite site = new BindingSite(String.valueOf(sequence, i, w), 
				//		i - mid, 
				//		scores[i], 
				//		strand);

				BindingSite site = new BindingSite(String.valueOf(sequence, i, w), 
						i, 
						scores[i], 
						strand);


				sites.add(site);

				//buffer = null;

				// skip ahead since we have covered this region
				i += w;
			}
		}

		return sites;
	}

	/*
	public static void processSite(char[] sequence,
			int start,
			int w,
			int mid,
			double score,
			char strand,
			List<BindingSite> sites) {
		//StringBuilder buffer = new StringBuilder();

		//for (int i = 0; i < w; ++i) {
		//	buffer.append(sequence[start + i]);
		//}

		Sequence seq = new Sequence(String.valueOf(sequence, start, w));

		//sites.add(new PotentialBindingSite(new Sequence(buffer.toString()), start - mid + 1, score));
		BindingSite site = new BindingSite(seq, 
				start - mid, 
				score, 
				strand);

		sites.add(site);
	}
	 */

	/**
	 * Return only the sites passing a score threshold.
	 * 
	 * @param sites
	 * @param threshold
	 * @return
	 */
	public static List<BindingSite> thresholdSites(List<BindingSite> sites, 
			double threshold) {
		List<BindingSite> ret = new ArrayList<BindingSite>();

		for (BindingSite site : sites) {
			//System.err.println(site.getScore() + " " + threshold + "\n");

			if (site.getScore() >= threshold) {
				ret.add(site);
			}
		}

		return ret;
	}

	/**
	 * Returns the binding site with the highest absolute score.
	 * 
	 * @param sites		The binding sites.
	 * @return			The binding site with the highest absolute score.
	 */
	public static BindingSite getBestBindingSites(List<BindingSite> sites) {
		double max = Double.MIN_VALUE;

		BindingSite best = null;

		for (BindingSite site : sites) {
			if (Math.abs(site.getScore()) > max) {
				max = Math.abs(site.getScore());
				best = site;
			}
		}

		return best;
	}

	/**
	 * Produce a version of the sequence where the binding sites
	 * are shown in upper case and the rest of the sequence in lower case.
	 * 
	 * @param sequence		The sequence.
	 * @param scores		The motif score at each base.
	 * @param threshold		The threshold to score a score significant.
	 * @return				The sequence with the binding sites in upper case.
	 */
	public static Sequence showScores(Sequence sequence, 
			double[] llkrf,
			double[] llkrr,
			double threshold) {
		String s = sequence.toString();

		StringBuilder buffer = new StringBuilder();

		for (int i = 0; i < sequence.getLength(); ++i) {
			if (llkrf[i] > threshold || llkrr[i] > threshold) {
				buffer.append(Character.toUpperCase(s.charAt(i)));
			} else {
				buffer.append(Character.toLowerCase(s.charAt(i)));
			}
		}

		return Sequence.create(sequence.getName(), buffer.toString());
	}

	public static Sequence showScores(char[] sequence, 
			BindingSite site,
			int n) {

		StringBuilder buffer = new StringBuilder();

		int mid = n / 2;

		int si = site.getOffset();
		int ei = si + site.getSequence().length();

		for (int i = 0; i < n; ++i) {
			int offset = i - mid;

			if (offset >= si && offset < ei) {
				buffer.append(Character.toUpperCase(sequence[i]));
			} else {
				buffer.append(Character.toLowerCase(sequence[i]));
			}
		}

		return Sequence.create(buffer.toString());
	}

	public static void search(Sequence sequence,
			Sequence revCompSequence,
			double[][] pwm,
			double bgscore,
			int w,
			final double[] llkrf,
			final double[] llkrr,
			double threshold) {

		char[] seq = sequence.toArray();
		char[] revCompSeq = revCompSequence.toArray();

		byte[] iSeq = Sequence.seqToIndexSeq(seq);
		byte[] iRevCompSeq = Sequence.seqToIndexSeq(revCompSeq);

		int n = seq.length;

		search(iSeq,
				iRevCompSeq,
				n,
				pwm,
				bgscore,
				w,
				threshold,
				llkrf,
				llkrr);
	}

	/**
	 * Search a sequence for motif binding sites. An array twice the
	 * length of the sequence is returned. The first n elements contain
	 * the scores on the forward strand, the second n elements contain
	 * the scores on the reverse strand.
	 * 
	 * @param sequence
	 * @param reverseComplement
	 * @param n
	 * @param motif
	 * @param bgscore
	 * @param w
	 * @param threshold
	 * @param llkrf
	 * @param llkrr
	 * @param fg
	 */
	public static void search(final byte[] sequence,
			final byte[] reverseComplement,
			int n,
			double[][] pwm,
			double bgscore,
			int w,
			double threshold,
			double[] llkrf,
			double[] llkrr) {

		search(sequence, n, pwm, bgscore, w, threshold, llkrf);

		search(reverseComplement, n, pwm, bgscore, w, threshold, llkrr);

		// Need to reverse llkrr so it is in the forward direction, since
		// the buffer may exceed the sequence length, only reverse the
		// first n bases.
		CollectionUtils.reverse(llkrr, 0, n);
	}

	/**
	 * Search for a motif and score all matches in the sequence.
	 * 
	 * @param sequence
	 * @param motif
	 * @param bgscore
	 * @param w
	 * @param scores
	 * @param threshold
	 */
	private static void search(byte[] sequence,
			int n,
			double[][] pwm,
			double bgscore,
			int w,
			double threshold,
			double[] scores) {
		//double[] fg = fglk(sequence, motif, w);
		fglk(sequence, n, pwm, w, FG);

		// the probability of an an equally weighted motif at this position
		score(n, FG, bgscore, threshold, scores);

		//return ret;
	}
	
	private static void search(byte[] sequence,
			Collection<Integer> startLocations,
			int n,
			double[][] pwm, 
			double bgscore,
			int w,
			double threshold,
			double[] scores) {
		//double[] fg = fglk(sequence, motif, w);
		fglk(sequence, startLocations, n, pwm, w, FG);

		// the probability of an an equally weighted motif at this position
		score(n, FG, bgscore, threshold, scores);

		//return ret;
	}

	/**
	 * Returns the theoretical max alignment score of a given motif.
	 * 
	 * @param motif
	 * @return
	 */
	public static double getMaxScore(Motif motif) {
		double p = 1;

		for (int i = 0; i < motif.getBaseCount(); ++i) {
			p *= motif.getCounts(i).getMaxScore();
		}

		double ratio = p / motif.getBgPwm();

		if (ratio > 1.0) {
			// Could use Math.log(ratio) to be like Homer
			return Mathematics.log2(ratio);
		} else {
			return 0;
		}
	}

	private static void score(int n,
			double[] fg, 
			double bgscore,
			double threshold,
			double[] scores) {
		double ratio;

		for (int i = 0; i < n; ++i) {
			ratio = fg[i] / bgscore; //bg[i];

			if (ratio > 1.0) {
				// Could use Math.log(ratio) to be like Homer
				scores[i] = Mathematics.log2(ratio);
			} else {
				scores[i] = 0;
			}

			// score must exceed a threshold to count
			if (scores[i] < threshold) {
				scores[i] = 0;
			}
		}

		//return ret;
	}

	/**
	 * 
	 * @param sequence
	 * @param n
	 * @param motif
	 * @param w				The motif width.
	 * @return
	 */
	private static void fglk(final byte[] sequence,
			int n,
			double[][] pwm,
			int w,
			double[] fg) {

		//double[] ret = Mathematics.zerosArray(n);

		// Need to reset fg to zero since it keeps a sum for the sequence
		CollectionUtils.fill(0, n, fg);
		//Arrays.fill(fg, 0, n, 0);

		int l = n - w + 1;

		byte base;

		// Move along the sequence
		for (int i = 0; i < l; ++i) {
			double p = 1;

			for (int j = 0; j < w; ++j) {
				//System.err.println(i + " " + j + " motif " + motif.getName() + " " + n + " " + w + " " + sequence.length);

				base = sequence[i + j];

				//System.err.print(base);

				p *= pwm[base][j]; //motif.getCount(base, j);

				//System.err.println(i + " " + j + " " + base + " " + Arrays.toString(pwm[base]));

				// No point continuing if p is zero, the matrix doesn't match
				// here and we learn nothing testing the other bases
				if (p == 0) {
					break;
				}
			}

			//System.err.println("p " + p);

			// Update all positions with the sum of the new p value if p > 0
			if (p > 0) {
				for (int j = 0; j < w; ++j) {
					fg[i + j] += p;
				}
			}
		}
	}

	private static void fglk(final byte[] sequence,
			Collection<Integer> startLocations,
			int n,
			double[][] pwm, 
			int w,
			double[] fg) {

		//double[] ret = Mathematics.zerosArray(n);

		// Need to reset fg to zero since it keeps a sum for the sequence
		CollectionUtils.fill(0, n, fg);
		//Arrays.fill(fg, 0, n, 0);

		byte base;
		
		int l = n - w;

		// skip around looking for triplets where we might match
		for (int i : startLocations) {
			// skip start locations that cannot contain the motif
			if (i > l) {
				continue;
			}
			
			double p = 1;

			for (int j = 0; j < w; ++j) {
				base = sequence[i + j];

				p *= pwm[base][j]; 

				// No point continuing if p is zero, the matrix doesn't match
				// here and we learn nothing testing the other bases
				if (p == 0) {
					break;
				}
			}

			if (p > 0) {
				//System.err.println("hmm " + p);
				for (int j = 0; j < w; ++j) {
					fg[i + j] += p;
				}
			}
		}
	}



	public static List<SearchRegion> getSearchRegions(Group group,
			int ext5p, 
			int ext3p,
			boolean mainVariants,
			boolean peakWidths,
			GenesDb genesDb) throws IOException, ParseException {
		List<SearchRegion> regions = new ArrayList<SearchRegion>();

		for (String id : group) {
			GenomicRegion region;

			if (GenomicRegion.isGenomicRegion(id)) {
				// assume it is a peak file that I made

				region = GenomicRegion.parse(id);

				if (peakWidths) {
					// Use the peak width as the search region rather than
					// extending around a midpoint
					regions.add(SearchRegion.createSearchRegion(region));
				} else {
					regions.add(SearchRegion.createSearchRegion(region, ext5p, ext3p));
				}
			} else {
				// See if the id is a symbol, in which case we are looking
				// at the TSS

				if (mainVariants) {
					regions.add(SearchRegion.createSearchRegion(genesDb.getMainGene(id), ext5p, ext3p));
				} else {
					List<Gene> genes = genesDb.getGenes(id);

					for (Gene gene : genes) {
						regions.add(SearchRegion.createSearchRegion(gene, ext5p, ext3p));
					}
				}
			}
		}

		return regions;
	}
}
