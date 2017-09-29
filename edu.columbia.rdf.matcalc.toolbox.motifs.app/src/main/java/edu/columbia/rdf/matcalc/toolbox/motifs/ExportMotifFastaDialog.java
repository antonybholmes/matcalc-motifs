package edu.columbia.rdf.matcalc.toolbox.motifs;

import java.awt.Dimension;

import javax.swing.Box;

import org.jebtk.modern.UI;
import org.jebtk.modern.button.ModernButtonGroup;
import org.jebtk.modern.button.ModernRadioButton;
import org.jebtk.modern.dialog.ModernDialogTaskWindow;
import org.jebtk.modern.panel.VBox;
import org.jebtk.modern.window.ModernWindow;
import org.jebtk.modern.window.WindowWidgetFocusEvents;

/**
 * Control which conservation scores are shown.
 * 
 * @author Antony Holmes Holmes
 *
 */
public class ExportMotifFastaDialog extends ModernDialogTaskWindow {
	private static final long serialVersionUID = 1L;

	private ModernRadioButton mCheckMotifs = 
			new ModernRadioButton("Motifs", true);
	
	private ModernRadioButton mCheckRegions = new ModernRadioButton("Regions");

	public ExportMotifFastaDialog(ModernWindow parent) {
		super(parent);

		setTitle("Export Motifs");

		setup();

		createUi();
	}

	private void setup() {
		addWindowListener(new WindowWidgetFocusEvents(mOkButton));

		setSize(new Dimension(320, 240));

		UI.centerWindowToScreen(this);
	}

	private final void createUi() {
		//this.getContentPane().add(new JLabel("Change " + getProductDetails().getProductName() + " settings", JLabel.LEFT), BorderLayout.PAGE_START);

		Box box = VBox.create();

		box.add(mCheckMotifs);

		box.add(UI.createVGap(10));

		box.add(mCheckRegions);

		setDialogCardContent(box);


		ModernButtonGroup group = new ModernButtonGroup();

		group.add(mCheckMotifs);
		group.add(mCheckRegions);
	}

	public boolean getMotifWise() {
		return mCheckMotifs.isSelected();
	}
}
