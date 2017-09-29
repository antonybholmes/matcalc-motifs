package edu.columbia.rdf.matcalc.toolbox.motifs.seqlogo;


import java.io.IOException;
import java.nio.file.Path;
import java.util.Collection;
import java.util.List;

import javax.xml.parsers.ParserConfigurationException;

import org.apache.batik.transcoder.TranscoderException;
import org.jebtk.bioinformatics.motifs.Motif;
import org.jebtk.bioinformatics.ui.BioInfDialog;
import org.jebtk.bioinformatics.ui.motifs.MotifModel;
import org.jebtk.bioinformatics.ui.motifs.MotifsTreePanel;
import org.jebtk.bioinformatics.ui.motifs.SeqLogosCanvas;
import org.jebtk.core.collections.CollectionUtils;
import org.jebtk.core.event.ChangeEvent;
import org.jebtk.core.event.ChangeListener;
import org.jebtk.core.io.FileUtils;
import org.jebtk.core.io.PathUtils;
import org.jebtk.graphplot.Image;
import org.jebtk.graphplot.MatrixGroupModel;
import org.jebtk.graphplot.ModernPlotCanvas;
import org.jebtk.graphplot.figure.FigurePanel;
import org.jebtk.modern.UI;
import org.jebtk.modern.UIService;
import org.jebtk.modern.contentpane.CenterTab;
import org.jebtk.modern.contentpane.CloseableHTab;
import org.jebtk.modern.contentpane.HTab;
import org.jebtk.modern.contentpane.ModernHContentPane;
import org.jebtk.modern.contentpane.SizableContentPane;
import org.jebtk.modern.dialog.DialogEvent;
import org.jebtk.modern.dialog.DialogEventListener;
import org.jebtk.modern.dialog.ModernDialogStatus;
import org.jebtk.modern.dialog.ModernMessageDialog;
import org.jebtk.modern.event.ModernClickEvent;
import org.jebtk.modern.event.ModernClickListener;
import org.jebtk.modern.event.ModernSelectionListener;
import org.jebtk.modern.graphics.colormap.ColorMapModel;
import org.jebtk.modern.graphics.icons.QuickOpenVectorIcon;
import org.jebtk.modern.graphics.icons.QuickSaveVectorIcon;
import org.jebtk.modern.help.ModernAboutDialog;
import org.jebtk.modern.help.RibbonPanelProductInfo;
import org.jebtk.modern.io.OpenRibbonPanel;
import org.jebtk.modern.io.RecentFilesService;
import org.jebtk.modern.io.SaveAsRibbonPanel;
import org.jebtk.modern.options.ModernOptionsRibbonPanel;
import org.jebtk.modern.panel.ModernPanel;
import org.jebtk.modern.ribbon.QuickAccessButton;
import org.jebtk.modern.ribbon.RibbonMenuItem;
import org.jebtk.modern.scrollpane.FixedIncScroller;
import org.jebtk.modern.scrollpane.ModernScrollPane;
import org.jebtk.modern.widget.ModernClickWidget;
import org.jebtk.modern.widget.ModernWidget;
import org.jebtk.modern.widget.tooltip.ModernToolTip;
import org.jebtk.modern.window.ModernRibbonWindow;
import org.jebtk.modern.zoom.ModernStatusZoomSlider;
import org.jebtk.modern.zoom.ZoomModel;
import org.jebtk.modern.zoom.ZoomRibbonSection;
import org.xml.sax.SAXException;


/**
 * Window for showing 2D graphs such as a scatter plot.
 * 
 * @author Antony Holmes Holmes
 *
 */
public class MainSeqLogoWindow extends ModernRibbonWindow implements ModernClickListener {
	private static final long serialVersionUID = 1L;

	private OpenRibbonPanel mOpenPanel = new OpenRibbonPanel();

	private SaveAsRibbonPanel mSaveAsPanel = new SaveAsRibbonPanel();

	protected ZoomModel mZoomModel = new ZoomModel();

	protected ColorMapModel mColorMapModel = 
			new ColorMapModel();

	protected MatrixGroupModel mGroupsModel = 
			new MatrixGroupModel();

	private SeqLogosCanvas mFigure = new SeqLogosCanvas(); //new PeaksCanvas(mBedGraphsModel, mGenomicModel);

	private MotifModel mMotifModel = new MotifModel();

	private FormatPanel mFormatPanel;

	private Path mFile;

	private MotifsTreePanel mMotifsPanel;

	//private ModernCheckBox mButtonRevComp;

	private MotifViewModel mViewModel = new MotifViewModel();

	private class MotifEvents implements ModernSelectionListener {

		@Override
		public void selectionChanged(ChangeEvent e) {
			openMotifs();
		}

	}

	private class MotifViewEvents implements ChangeListener {

		@Override
		public void changed(ChangeEvent e) {
			openMotifs();
		}
	}

	private class ExportCallBack implements DialogEventListener {
		private Path mFile;
		private Path mPwd;

		public ExportCallBack(Path file, Path pwd) {
			mFile = file;
			mPwd = pwd;
		}

		@Override
		public void statusChanged(DialogEvent e) {
			if (e.getStatus() == ModernDialogStatus.OK) {
				try {	
					save(mFile);
				} catch (Exception ex) {
					ex.printStackTrace();
				}
			} else {
				try {
					export(mPwd);
				} catch (Exception ex) {
					ex.printStackTrace();
				}
			}
		}
	}

	public MainSeqLogoWindow() {
		super(new SeqLogoInfo());
		
		init();
	}

	public MainSeqLogoWindow(Path file) {
		this();

		try {
			openFile(file);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public MainSeqLogoWindow(Collection<String> names) {
		this();

		mMotifsPanel.search(names);
	}

	public void init() {
		mMotifsPanel = new MotifsTreePanel(this, mMotifModel);

		createRibbon();

		createUi();

		mMotifModel.addSelectionListener(new MotifEvents());

		//mButtonRevComp.addClickListener(new RevCompEvents());

		mViewModel.addChangeListener(new MotifViewEvents());

		setSize(1280, 720);

		UI.centerWindowToScreen(this);
	}

	public final void createRibbon() {
		//RibbongetRibbonMenu() getRibbonMenu() = new RibbongetRibbonMenu()(0);
		RibbonMenuItem menuItem;

		menuItem = new RibbonMenuItem(UI.MENU_OPEN);
		getRibbonMenu().addTabbedMenuItem(menuItem, mOpenPanel);

		menuItem = new RibbonMenuItem(UI.MENU_SAVE_AS);
		getRibbonMenu().addTabbedMenuItem(menuItem, mSaveAsPanel);

		menuItem = new RibbonMenuItem(UI.MENU_EXIT);
		getRibbonMenu().addTabbedMenuItem(menuItem);

		getRibbonMenu().addSeparator();

		menuItem = new RibbonMenuItem(UI.MENU_INFO);
		getRibbonMenu().addTabbedMenuItem(menuItem, 
				new RibbonPanelProductInfo(getAppInfo()));
		//getRibbonMenu().addTabbedMenuItem(menuItem);

		menuItem = new RibbonMenuItem(UI.MENU_OPTIONS);
		getRibbonMenu().addTabbedMenuItem(menuItem, new ModernOptionsRibbonPanel(getAppInfo()));


		getRibbonMenu().setDefaultIndex(0);

		getRibbonMenu().addClickListener(this);


		ModernClickWidget button;

		//Ribbon2 ribbon = new Ribbon2();
		getRibbon().setHelpButtonEnabled(getAppInfo());

		button = new QuickAccessButton(UIService.getInstance().loadIcon(QuickOpenVectorIcon.class, 16));
		button.setClickMessage(UI.MENU_OPEN);
		button.setToolTip(new ModernToolTip("Open", 
				"Open peak files."));
		button.addClickListener(this);
		addQuickAccessButton(button);

		button = new QuickAccessButton(UIService.getInstance().loadIcon(QuickSaveVectorIcon.class, 16));
		button.setClickMessage(UI.MENU_SAVE);
		button.setToolTip(new ModernToolTip("Save", 
				"Save the current image."));
		button.addClickListener(this);
		addQuickAccessButton(button);



		// home
		//RibbonToolbar toolbar = new RibbonToolbar(getRibbon(), "Plot");

		//RibbonSection toolbarContainer = new RibbonSection("Options");

		//mButtonRevComp = new RibbonLargeCheckButtonIconText2("Reverse Complement", 
		//		UIResources.getInstance().loadIcon("reverse_complement", 24),
		//		"Reverse Complement",
		//		"Reverse complement motifs.");

		//Ui.setSize(mButtonRevComp, Ribbon2.MEDIUM_TEXT_BUTTON_SIZE);

		getRibbon().getToolbar("Plot").add(new MotifViewRibbonSection(getRibbon(), mViewModel));

		//mButtonRevComp = new ModernCheckBox("Reverse Complement");
		//toolbarContainer.add(new RibbonStripContainer(mButtonRevComp));
		//toolbar.add(toolbarContainer);
		
		BaseColorRibbonSection baseColorRibbon = new BaseColorRibbonSection(this);
		
		baseColorRibbon.addClickListener(new ModernClickListener(){

			@Override
			public void clicked(ModernClickEvent e) {
				openMotifs();
			}});
		
		getRibbon().getToolbar("Plot").add(baseColorRibbon);

		//toolbarSection = new PlotSizeRibbonSection(mCanvas.getGraphSpace().getLayoutProperties());
		//toolbar.add(toolbarSection);

		//ZoomRibbonSection zoomSection = 
		//		new ZoomRibbonSection(this, zoomModel, ribbon);

		//toolbar.add(zoomSection);

		//LegendRibbonSection legendSection =
		//		new LegendRibbonSection(mCanvas.getGraphProperties().getLegend());

		//toolbar.add(legendSection);



		//getRibbon().addToolbar(toolbar);

		getRibbon().getToolbar("View").add(new ZoomRibbonSection(this, mZoomModel));

		//setRibbon(ribbon, getRibbonMenu());

		getRibbon().setSelectedIndex(1);
	}

	public final void createUi() {

		//ImageCanvas imageCanvas = new ImageCanvas(mCanvas);

		//ZoomCanvas zoomCanvas = new ZoomCanvas(canvas);

		ModernPlotCanvas canvas = new FigurePanel(mFigure);
		
		canvas.setZoomModel(mZoomModel);

		ModernScrollPane scrollPane = new ModernScrollPane(canvas);
		
		scrollPane.getVScrollBar().setScroller(new FixedIncScroller(200));

		ModernPanel panel = new ModernPanel(scrollPane, ModernWidget.BORDER);

		setCard(panel);

		mStatusBar.addRight(new ModernStatusZoomSlider(mZoomModel));
		
		addMotifPanel();
	}

	private void addMotifPanel() {
		if (getTabsPane().getModel().getLeftTabs().containsTab("Motifs")) {
			return;
		}

		mMotifsPanel.setBorder(ModernPanel.DOUBLE_BORDER);

		getTabsPane().addLeftTab("Motifs", 
				new HTab("Motifs", mMotifsPanel), 250, 200, 500);
	}

	private void addFormatPane() {
		if (getTabsPane().getModel().getRightTabs().containsTab("Format Plot")) {
			return;
		}

		getTabsPane().addRightTab("Format Plot", 
				new CloseableHTab("Format Plot", mFormatPanel, getTabsPane()), 300, 200, 500);

	}

	@Override
	public final void clicked(ModernClickEvent e) {
		if (e.getMessage().equals(UI.MENU_OPEN) ||
				e.getMessage().equals(UI.MENU_BROWSE) ||
				e.getMessage().startsWith("Other...")) {
			try {
				browseForFile();
			} catch (IOException e1) {
				e1.printStackTrace();
			} catch (SAXException e1) {
				e1.printStackTrace();
			} catch (ParserConfigurationException e1) {
				e1.printStackTrace();
			}
		} else if (e.getMessage().equals(OpenRibbonPanel.FILE_SELECTED)) {
			try {
				openFile(mOpenPanel.getSelectedFile());
			} catch (IOException e1) {
				e1.printStackTrace();
			}
		} else if (e.getMessage().equals(OpenRibbonPanel.DIRECTORY_SELECTED)) {
			try {
				browseForFile(mOpenPanel.getSelectedDirectory());
			} catch (IOException e1) {
				e1.printStackTrace();
			} catch (SAXException e1) {
				e1.printStackTrace();
			} catch (ParserConfigurationException e1) {
				e1.printStackTrace();
			}
		} else if (e.getMessage().equals(UI.MENU_SAVE)) {
			try {
				export();
			} catch (IOException e1) {
				e1.printStackTrace();
			} catch (TranscoderException e1) {
				e1.printStackTrace();
			}
		} else if (e.getMessage().equals(SaveAsRibbonPanel.DIRECTORY_SELECTED)) {
			try {
				export(mSaveAsPanel.getSelectedDirectory());
			} catch (IOException e1) {
				e1.printStackTrace();
			} catch (TranscoderException e1) {
				e1.printStackTrace();
			}
		} else if (e.getMessage().equals("Format Plot")) {
			addFormatPane();
		} else if (e.getMessage().equals(UI.MENU_ABOUT)) {
			ModernAboutDialog.show(this, getAppInfo());
		} else if (e.getMessage().equals(UI.MENU_EXIT)) {
			close();
		} else {

		}
	}

	private void browseForFile() throws IOException, SAXException, ParserConfigurationException {
		browseForFile(RecentFilesService.getInstance().getPwd());
	}

	private void browseForFile(Path workingDirectory) throws IOException, SAXException, ParserConfigurationException {
		openFile(BioInfDialog.openMotifFile(this, workingDirectory));
	}

	public void openFile(Path file) throws IOException {
		if (file == null) {
			return;
		}

		if (mFile != null) {
			// We have already opened a file so create a new
			// window

			MainSeqLogoWindow window = new MainSeqLogoWindow(file);

			window.setVisible(true);
		} else {

			if (PathUtils.getFileExt(file).equals("motif")) {
				Motif motif = Motif.parseMotif(file);

				openMotif(motif);
			} else if (PathUtils.getFileExt(file).equals("pwm")) {
				Motif motif = Motif.parsePwmMotif(file);

				openMotif(motif);
			} else {
				// do nothing
			}

			RecentFilesService.getInstance().add(file);

			mFile = file;
		}
	}

	private void openMotifs() {
		openMotifs(mMotifModel.getItems());
	}

	private void openMotif(Motif motif) {
		openMotifs(CollectionUtils.asList(motif));
	}

	private void openMotifs(List<Motif> motifs) {
		if (mViewModel.getRevComp()) {
			mFigure.setMotifs(Motif.reverseComplement(motifs), mViewModel.get());
		} else {
			mFigure.setMotifs(motifs, mViewModel.get());
		}

		//System.err.println("open motif " + motifs.get(0).getName());
	}


	private void export() throws IOException, TranscoderException {
		export(RecentFilesService.getInstance().getPwd());
	}

	private void export(Path pwd) throws IOException, TranscoderException {
		Path file = Image.saveFile(this, pwd);

		if (file == null) {
			return;
		}

		if (FileUtils.exists(file)) {
			ModernMessageDialog.createFileReplaceDialog(this, file, new ExportCallBack(file, pwd));
		} else {
			save(file);
		}
	}

	private void save(Path file) throws IOException, TranscoderException {
		Image.write(mFigure, file);

		ModernMessageDialog.createFileSavedDialog(this, file);
	}
}
