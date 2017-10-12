package edu.columbia.rdf.matcalc.toolbox.motifs;

import java.io.BufferedWriter;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import javax.swing.JFrame;

import org.jebtk.bioinformatics.dna.GenomeAssemblyLocal;
import org.jebtk.bioinformatics.dna.GenomeAssemblyWeb;
import org.jebtk.bioinformatics.dna.Sequence;
import org.jebtk.bioinformatics.genomic.Chromosome;
import org.jebtk.bioinformatics.genomic.ChromosomeSizesService;
import org.jebtk.bioinformatics.genomic.GenomeAssembly;
import org.jebtk.bioinformatics.genomic.GenomeAssemblyService;
import org.jebtk.bioinformatics.genomic.GenomicRegion;
import org.jebtk.bioinformatics.motifs.Motif;
import org.jebtk.bioinformatics.motifs.MotifsDataSourceService;
import org.jebtk.bioinformatics.motifs.MotifsFs;
import org.jebtk.bioinformatics.motifs.MotifsXmlFs;
import org.jebtk.bioinformatics.ui.external.ucsc.BedGuiFileFilter;
import org.jebtk.bioinformatics.ui.filters.MotifGuiFileFilter;
import org.jebtk.bioinformatics.ui.filters.MotifPwmGuiFileFilter;
import org.jebtk.bioinformatics.ui.groups.GroupsModel;
import org.jebtk.bioinformatics.ui.groups.GroupsPanel;
import org.jebtk.core.Resources;
import org.jebtk.core.collections.CollectionUtils;
import org.jebtk.core.io.FileUtils;
import org.jebtk.core.io.PathUtils;
import org.jebtk.core.settings.SettingsService;
import org.jebtk.core.text.Join;
import org.jebtk.math.matrix.DataFrame;
import org.jebtk.modern.UI;
import org.jebtk.modern.UIService;
import org.jebtk.modern.contentpane.CloseableHTab;
import org.jebtk.modern.contentpane.SizableContentPane;
import org.jebtk.modern.dialog.MessageDialogType;
import org.jebtk.modern.dialog.ModernDialogStatus;
import org.jebtk.modern.dialog.ModernMessageDialog;
import org.jebtk.modern.event.ModernClickEvent;
import org.jebtk.modern.event.ModernClickListener;
import org.jebtk.modern.graphics.icons.SaveVectorIcon;
import org.jebtk.modern.graphics.icons.SearchVectorIcon;
import org.jebtk.modern.help.GuiAppInfo;
import org.jebtk.modern.io.FileDialog;
import org.jebtk.modern.io.GuiFileExtFilter;
import org.jebtk.modern.io.RecentFilesService;
import org.jebtk.modern.ribbon.Ribbon;
import org.jebtk.modern.ribbon.RibbonLargeButton;
import org.jebtk.modern.widget.ModernClickWidget;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import edu.columbia.rdf.matcalc.FileType;
import edu.columbia.rdf.matcalc.MainMatCalcWindow;
import edu.columbia.rdf.matcalc.toolbox.CalcModule;
import edu.columbia.rdf.matcalc.toolbox.motifs.app.MotifsInfo;
import edu.columbia.rdf.matcalc.toolbox.motifs.seqlogo.MainSeqLogoWindow;
import edu.columbia.rdf.matcalc.toolbox.motifs.seqlogo.SeqLogoIcon;


public class MotifsModule extends CalcModule implements ModernClickListener {
	public static final Logger LOG = 
			LoggerFactory.getLogger(MotifsModule.class);

	private static final Path CHR_SIZE_FILE = 
			PathUtils.getPath("res/modules/motifs/ucsc_chromosome_sizes_hg19.txt.gz");

	private static final Path DATABASE_DIR = 
			PathUtils.getPath("res/modules/motifs/database");
	
	private static final Path RES_DIR = PathUtils.getPath("res/modules/dna");
	
	
	private static final Set<GuiFileExtFilter> FILE_TYPES_SET =
			new TreeSet<GuiFileExtFilter>();

	private MainMatCalcWindow mWindow;

	private GroupsModel mRegionGroupsModel = new GroupsModel();
	private GroupsPanel mRegionGroupsPanel;


	//private GenomeAssemblyWeb mAssembly;

	static {
		if (SettingsService.getInstance().getAsBool("org.matcalc.toolbox.bio.dna.web.enabled")) {
			try {
				GenomeAssemblyService.getInstance().add(new GenomeAssemblyWeb(new URL(SettingsService.getInstance().getAsString("dna.remote-url"))));
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		
		GenomeAssemblyService.getInstance().add(new GenomeAssemblyLocal(RES_DIR));

		FILE_TYPES_SET.add(new MotifGuiFileFilter());
		FILE_TYPES_SET.add(new MotifPwmGuiFileFilter());
	}

	public MotifsModule() throws MalformedURLException, IOException {
		//GenomeAssemblyService.getInstance().add(new GenomeAssemblyWeb(new URL(SettingsService.getInstance().getAsString("dna.remote-url"))));
		//GenomeAssemblyService.getInstance().add(new GenomeAssemblyFile4Bit(FILE_PATH));

		//mAssembly = new GenomeAssemblyWeb(SettingsService.getInstance().getAsUrl("dna.remote-url"));

		//MotifsDBService.getInstance().addBackEnd(new MotifsWeb());

		MotifsDataSourceService.getInstance().addDataSource(new MotifsFs(DATABASE_DIR));
		MotifsDataSourceService.getInstance().addDataSource(new MotifsXmlFs(DATABASE_DIR));
	}

	@Override
	public void run(String... args) {
		// Do nothing
	}

	@Override
	public String getName() {
		return "Motifs";
	}

	@Override
	public GuiAppInfo getModuleInfo() {
		return new MotifsInfo();
	}

	@Override
	public Set<GuiFileExtFilter> getOpenFileFilters() {
		return FILE_TYPES_SET;
	}

	@Override
	public void init(MainMatCalcWindow window) {
		mWindow = window;

		mRegionGroupsPanel = new GroupsPanel(mWindow, mRegionGroupsModel);

		Ribbon ribbon = window.getRibbon();

		ModernClickWidget button;

		/*
		button = new RibbonLargeButton2("Motif Groups", 
				UIResources.getInstance().loadIcon("sidebar", 24),
				"Motif Groups", 
				"Show the motifs group");
		button.addClickListener(this);
		ribbon.getToolbar("DNA").getSection("Motifs").add(button);
		 */

		button = new RibbonLargeButton("Search", 
				UIService.getInstance().loadIcon(SearchVectorIcon.class, 24),
				"Search", 
				"Search for motifs.");
		button.addClickListener(this);
		ribbon.getToolbar("DNA").getSection("Motifs").add(button);

		ribbon.getToolbar("DNA").getSection("Motifs").addSeparator();

		button = new RibbonLargeButton("GC Background", 
				UIService.getInstance().loadIcon("enrichment", 24),
				"GC Background", 
				"Generate background.");
		button.addClickListener(this);
		ribbon.getToolbar("DNA").getSection("Motifs").add(button);

		button = new RibbonLargeButton("Enrichment", 
				UIService.getInstance().loadIcon("enrichment", 24),
				"Enrichment", 
				"Look for enriched motifs.");
		button.addClickListener(this);
		ribbon.getToolbar("DNA").getSection("Motifs").add(button);

		ribbon.getToolbar("DNA").getSection("Motifs").addSeparator();

		button = new RibbonLargeButton("Export BED", 
				UIService.getInstance().loadIcon(SaveVectorIcon.class, 24),
				"Export BED", 
				"Export Results as BED.");
		button.addClickListener(this);
		ribbon.getToolbar("DNA").getSection("Motifs").add(button);

		ribbon.getToolbar("DNA").getSection("Motifs").addSeparator();

		button = new RibbonLargeButton("SeqLogo", 
				UIService.getInstance().loadIcon(SeqLogoIcon.class, 24),
				"SeqLogo", 
				"Browse sequence logos.");
		button.addClickListener(this);
		ribbon.getToolbar("DNA").getSection("Motifs").add(button);
	}

	@Override
	public void clicked(ModernClickEvent e) {
		if (e.getMessage().equals("Export BED")) {
			try {
				exportBed();
			} catch (IOException e1) {
				e1.printStackTrace();
			}
		} else if (e.getMessage().equals("Motif Groups")) {
			addRegionGroupsPanel();
		} else if (e.getMessage().equals("Search")) {
			try {
				motifSearch();
			} catch (IOException e1) {
				e1.printStackTrace();
			}
		} else if (e.getMessage().equals("Enrichment")) {
			try {
				motifEnrichment();
			} catch (IOException e1) {
				e1.printStackTrace();
			}
		} else if (e.getMessage().equals("GC Background")) {
			try {
				gcBackground();
			} catch (IOException e1) {
				e1.printStackTrace();
			}
		} else if (e.getMessage().equals("SeqLogo")) {
			try {
				seqLogo();
			} catch (IOException e1) {
				e1.printStackTrace();
			}
		} else {
			// do nothing
		}
	}

	@Override
	public DataFrame autoOpenFile(final MainMatCalcWindow window,
			final Path file,
			FileType type, 
			int headers,
			int rowAnnotations,
			String delimiter,
			Collection<String> skipLines) throws IOException {
		seqLogo(file);

		return null;
	}

	private void exportBed() throws IOException {
		Path file = FileDialog.save(mWindow).filter(new BedGuiFileFilter()).getFile(RecentFilesService.getInstance().getPwd());

		if (file == null) {
			return;
		}

		if (FileUtils.exists(file)) {
			if (ModernMessageDialog.createFileReplaceDialog(mWindow, file) == ModernDialogStatus.CANCEL) {
				return;
			}
		}

		DataFrame m = mWindow.getCurrentMatrix();

		int locCol = DataFrame.findColumn(m, "Motif Location");
		int scoreCol = DataFrame.findColumn(m, "Score");
		int strandCol = DataFrame.findColumn(m, "Strand");
		int nameCol = DataFrame.findColumn(m, "Motif Name");
		int idCol = DataFrame.findColumn(m, "Motif ID");

		BufferedWriter writer = FileUtils.newBufferedWriter(file);

		try {
			writer.write("track name=\"Motifs\" description=\"Motifs\"");
			writer.newLine();

			for (int i = 0; i < m.getRowCount(); ++i) {
				GenomicRegion r = GenomicRegion.parse(m.getText(i, locCol));
				double score = m.getValue(i, scoreCol);
				char strand = m.getText(i, strandCol).charAt(0);

				String name = m.getText(i, nameCol).replaceAll(" +", "_").replaceAll(",", "") + 
						"_" + 
						m.getText(i, idCol);

				writer.write(Join.onTab().values(r.getChr(), 
						r.getStart(), 
						r.getEnd(),
						name,
						score,
						strand).toString());
				writer.newLine();
			}
		} finally {
			writer.close();
		}

		ModernMessageDialog.createFileSavedDialog(mWindow, file);
	}

	private void seqLogo() throws IOException {
		// Lets determine if we have a motif table loaded, and if so can
		// we extract some motif names and default to showing those first.

		int c = -1;
		
		DataFrame m = mWindow.getCurrentMatrix();

		if (m != null) {
			c = DataFrame.findColumn(m, "Motif ID");
		}

		JFrame window;

		if (c != -1) {
			List<String> names = CollectionUtils.uniquePreserveOrder(m.columnAsText(c));

			window = new MainSeqLogoWindow(names);

		} else {
			window = new MainSeqLogoWindow();
		}

		UI.centerWindowToScreen(window);

		window.setVisible(true);
	}

	private void seqLogo(Path file) throws IOException {
		JFrame window = new MainSeqLogoWindow(file);

		UI.centerWindowToScreen(window);

		window.setVisible(true);
	}

	private void addRegionGroupsPanel() {
		if (mWindow.getTabsPane().getModel().getLeftTabs().containsTab("Motif Groups")) {
			return;
		}

		SizableContentPane sizePane = new SizableContentPane("Motif Groups", 
				new CloseableHTab("Motif Groups", mRegionGroupsPanel, mWindow.getTabsPane()), 
				250, 
				100, 
				500);

		mWindow.getTabsPane().getModel().addLeftTab(sizePane);
	}

	private void motifSearch() throws IOException {
		MotifSearchDialog dialog = 
				new MotifSearchDialog(mWindow, mRegionGroupsModel);

		dialog.setVisible(true);

		if (dialog.getStatus() == ModernDialogStatus.CANCEL) {
			return;
		}

		MotifSearchTask task = new MotifSearchTask(mWindow,
				dialog.getMotifs(),
				dialog.getThreshold());

		task.doInBackground();
		task.done();
	}

	private void motifEnrichment() throws IOException {
		MotifEnrichmentDialog dialog = new MotifEnrichmentDialog(mWindow);

		dialog.setVisible(true);

		if (dialog.getStatus() == ModernDialogStatus.CANCEL) {
			return;
		}

		ModernDialogStatus status = ModernMessageDialog.createDialog(mWindow,
				"Searching for motifs may take several minutes.",
				MessageDialogType.INFORMATION_OK_CANCEL);

		if (status == ModernDialogStatus.CANCEL) {
			return;
		}

		List<Motif> searchMotifs = dialog.getMotifs();

		double threshold = dialog.getThreshold();
		double sensitivity = dialog.getMinSensitivity();
		double specificity = dialog.getMinSpecificity();
		MainMatCalcWindow foregroundGroup = mWindow; //dialog.getForegroundGroup();
		MainMatCalcWindow backgroundGroup = dialog.getBackgroundGroup();

		// Load some chromosome sizes for hg19

		ChromosomeSizesService.getInstance().load(GenomeAssembly.HG19, 
				Resources.getGzipReader(CHR_SIZE_FILE));

		if (dialog.useForeVsBackMode()) {
			MotifEnrichmentTask task = new MotifEnrichmentTask(mWindow,
					searchMotifs,
					foregroundGroup,
					backgroundGroup,
					threshold,
					sensitivity,
					specificity);

			task.doInBackground(); //execute();
		} else {
			MotifEnrichmentGCHistTask task = new MotifEnrichmentGCHistTask(mWindow,
					"hg19",
					GenomeAssemblyService.getInstance(),
					ChromosomeSizesService.getInstance().getSizes(GenomeAssembly.HG19),
					searchMotifs,
					foregroundGroup,
					threshold,
					sensitivity,
					specificity);

			task.doInBackground();
		}
	}

	private void gcBackground() throws IOException {
		ModernDialogStatus status = ModernMessageDialog.createDialog(mWindow,
				"Generating a GC matched set of sequences may take several minutes.",
				MessageDialogType.INFORMATION_OK_CANCEL);

		if (status == ModernDialogStatus.CANCEL) {
			return;
		}

		// Load some chromosome sizes for hg19
		ChromosomeSizesService.getInstance().load(GenomeAssembly.HG19, 
				Resources.getGzipReader(CHR_SIZE_FILE));

		GCBackgroundTask task = new GCBackgroundTask(mWindow,
				GenomeAssembly.HG19,
				GenomeAssemblyService.getInstance(),
				ChromosomeSizesService.getInstance().getSizes(GenomeAssembly.HG19));

		task.doInBackground();
	}

	public static List<SearchSequence> matrixToSequences(DataFrame m) {
		int dnaLocationColumn = -1;
		int chrCol = -1;
		int startCol = -1;
		int endCol = -1;
		int idCol = -1;

		// Find a column refering to a genomic location
		for (int i = 0; i < m.getColumnCount(); ++i) {
			if (GenomicRegion.isGenomicRegion(m.getText(0, i))) {
				dnaLocationColumn = i;
				break;
			}
		}

		// If this is not found test whether we have chr, start and end columns
		if (dnaLocationColumn == -1) {
			for (int i = 0; i < m.getColumnCount(); ++i) {
				if (Chromosome.isChr(m.getText(0, i))) {
					chrCol = i;
					break;
				}
			}

			if (chrCol != -1) {
				startCol = chrCol + 1;
				endCol = chrCol + 2;
			}
		}
		
		// Last resort pick the first column as an id
		if (chrCol == -1) {
			idCol = 0;
		}

		int dnaColumn = -1;

		for (int i = 0; i < m.getColumnCount(); ++i) {
			if (isDna(m.getText(0, i))) {
				dnaColumn = i;
				break;
			}
		}

		if (dnaColumn == -1) {
			return Collections.emptyList();
		}

		List<SearchSequence> sequences = 
				new ArrayList<SearchSequence>(m.getRowCount());

		if (dnaLocationColumn != -1) {
			for (int i = 0; i < m.getRowCount(); ++i) {
				GenomicRegion region = GenomicRegion.parse(m.getText(i, dnaLocationColumn));
				String dna = m.getText(i, dnaColumn);
				
				sequences.add(new SearchSequence(region, Sequence.create(dna)));
			}
		} else if (dnaLocationColumn != -1){
			for (int i = 0; i < m.getRowCount(); ++i) {
				GenomicRegion region = GenomicRegion.create(Chromosome.parse(m.getText(i, chrCol)), 
						(int)m.getValue(i, startCol), 
						(int)m.getValue(i, endCol));
				
				String dna = m.getText(i, dnaColumn);
				
				sequences.add(new SearchSequence(region, Sequence.create(dna)));
			}
		} else {
			for (int i = 0; i < m.getRowCount(); ++i) {
				String id = m.getText(i, idCol);
				
				String dna = m.getText(i, dnaColumn);
				
				sequences.add(new SearchSequence(id, Sequence.create(dna)));
			}
		}

		return sequences;
	}

	/**
	 * The DNA sequence must be at least 10 bp long to be considered useful.
	 * This is to stop short labels in columns such as 'a' from being
	 * misinterpreted as DNA.
	 * 
	 * @param text
	 * @return
	 */
	private static boolean isDna(String text) {
		return Sequence.DNA_REGEX.matcher(text).matches() && text.length() > 10;
	}


}
