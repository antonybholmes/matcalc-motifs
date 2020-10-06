package edu.columbia.rdf.matcalc.toolbox.motifs;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.jebtk.bioinformatics.genomic.Gene;
import org.jebtk.bioinformatics.genomic.Genome;
import org.jebtk.bioinformatics.genomic.GenomicRegion;
import org.jebtk.bioinformatics.genomic.SequenceReader;
import org.jebtk.bioinformatics.genomic.SequenceRegion;
import org.jebtk.bioinformatics.genomic.Strand;

/**
 * Describes what we looked for
 * 
 * @author Antony Holmes
 *
 */
public class SearchRegion {

  private GenomicRegion mRegion;
  private Strand mStrand;
  private GenomicRegion mReferenceRegion;
  private String mName;
  private int mExt5p;
  private int mExt3p;
  private GenomicRegion mSearchRegion;

  public SearchRegion(String name, GenomicRegion region, GenomicRegion referenceRegion, Strand strand, int ext5p,
      int ext3p) {
    mName = name;
    mRegion = region;
    mReferenceRegion = referenceRegion;

    mStrand = strand;
    mExt5p = ext5p;
    mExt3p = ext3p;

    if (strand == Strand.SENSE) {
      mSearchRegion = GenomicRegion.extend(referenceRegion, ext5p, ext3p);
    } else {
      mSearchRegion = GenomicRegion.extend(referenceRegion, ext3p, ext5p);
    }
  }

  public String getName() {
    return mName;
  }

  public int getExt5p() {
    return mExt5p;
  }

  public int getExt3p() {
    return mExt3p;
  }

  /**
   * The point around which searches are conducted, e.g. the mid point of a peak
   * or the TSS of a gene.
   * 
   * @return
   */
  public GenomicRegion getReferencePoint() {
    return mReferenceRegion;
  }

  // public double getScore() {
  // return mScore;
  // }

  public GenomicRegion getRegion() {
    return mRegion;
  }

  /**
   * The bounds of the search region 5' to 3' on the forward strand.
   * 
   * @return
   */
  public GenomicRegion getSearchRegion() {
    return mSearchRegion;
  }

  public Strand getStrand() {
    return mStrand;
  }

  public static List<SequenceRegion> getSequences(Genome genome, SequenceReader reader,
      List<SearchRegion> searchRegions) throws IOException {
    List<SequenceRegion> sequences = new ArrayList<SequenceRegion>();

    for (SearchRegion searchRegion : searchRegions) {
      sequences.add(reader.getSequence(genome, searchRegion.mSearchRegion));
    }

    return sequences;
  }

  /**
   * Create a peak search region using the peak width
   * 
   * @param region
   */
  public static SearchRegion createSearchRegion(GenomicRegion region) {
    GenomicRegion mid = GenomicRegion.midRegion(region);

    int ext5p = mid.getStart() - region.getStart();
    int ext3p = region.getEnd() - mid.getStart();

    return createSearchRegion(region, ext5p, ext3p);
  }

  /**
   * Search for a peak.
   * 
   * @param region
   */
  public static SearchRegion createSearchRegion(GenomicRegion region, int ext5p, int ext3p) {
    return new SearchRegion(region.getLocation(), region, GenomicRegion.midRegion(region), Strand.SENSE, ext5p, ext3p);
  }

  /**
   * Search for a gene.
   * 
   * @param gene
   */
  public static SearchRegion createSearchRegion(Gene gene, int ext5p, int ext3p) {
    return new SearchRegion(gene.getSymbol() + " (" + gene.getRefSeq() + ")", gene.getTss(), gene.getTss(),
        gene.getStrand(), ext5p, ext3p);
  }
}
