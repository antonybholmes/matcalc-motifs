package edu.columbia.rdf.matcalc.toolbox.motifs.seqlogo;

import org.jebtk.bioinformatics.genomic.DnaService;
import org.jebtk.modern.UI;
import org.jebtk.modern.button.ModernButton;
import org.jebtk.modern.dialog.ModernDialogFlatButton;
import org.jebtk.modern.dialog.ModernDialogHelpWindow;
import org.jebtk.modern.dialog.ModernDialogStatus;
import org.jebtk.modern.dialog.ModernMessageDialog;
import org.jebtk.modern.event.ModernClickEvent;
import org.jebtk.modern.event.ModernClickListener;
import org.jebtk.modern.graphics.color.ColorSwatchButton;
import org.jebtk.modern.panel.HExpandBox;
import org.jebtk.modern.panel.VBox;
import org.jebtk.modern.window.ModernWindow;
import org.jebtk.modern.window.WindowWidgetFocusEvents;

/**
 * Options for motif searching.
 * 
 * @author Antony Holmes Holmes
 *
 */
public class BaseDialog extends ModernDialogHelpWindow {
  private static final long serialVersionUID = 1L;

  private ColorSwatchButton mColorAButton;

  private ColorSwatchButton mColorCButton;

  private ColorSwatchButton mColorGButton;

  private ColorSwatchButton mColorTButton;

  private ColorSwatchButton mColorNButton;

  public BaseDialog(ModernWindow parent) {
    super(parent, "org.matcalc.toolbox.bio.motifs.seqlogo.dna-bases.help.url");

    setTitle("DNA Bases");

    setup();

    createUi();
  }

  private void setup() {
    addWindowListener(new WindowWidgetFocusEvents(mOkButton));

    setSize(400, 380);

    UI.centerWindowToScreen(this);
  }

  private final void createUi() {
    mColorAButton = new ColorSwatchButton(mParent,
        DnaService.getInstance().getBaseAColor());

    mColorCButton = new ColorSwatchButton(mParent,
        DnaService.getInstance().getBaseCColor());

    mColorGButton = new ColorSwatchButton(mParent,
        DnaService.getInstance().getBaseGColor());

    mColorTButton = new ColorSwatchButton(mParent,
        DnaService.getInstance().getBaseTColor());

    mColorNButton = new ColorSwatchButton(mParent,
        DnaService.getInstance().getBaseNColor());

    VBox box = VBox.create();

    sectionHeader("Base Color", box);

    box.add(new HExpandBox("A", mColorAButton));
    box.add(UI.createVGap(5));
    box.add(new HExpandBox("C", mColorCButton));
    box.add(UI.createVGap(5));
    box.add(new HExpandBox("G", mColorGButton));
    box.add(UI.createVGap(5));
    box.add(new HExpandBox("T", mColorTButton));
    box.add(UI.createVGap(5));
    box.add(new HExpandBox("N", mColorNButton));
    box.add(UI.createVGap(40));

    ModernButton button = new ModernDialogFlatButton("Defaults");
    button.addClickListener(new ModernClickListener() {

      @Override
      public void clicked(ModernClickEvent e) {
        resetToDefaults();
      }
    });

    box.add(button);

    setCard(box);
  }

  /**
   * Reset the colors to their defaults.
   */
  private void resetToDefaults() {
    ModernDialogStatus status = ModernMessageDialog.createOkCancelWarningDialog(
        mParent,
        "The base colors will be reset to their default values.");

    if (status == ModernDialogStatus.OK) {
      DnaService.getInstance().reset();

      setColors();
    }
  }

  /**
   * Change the button color.
   */
  private void setColors() {
    mColorAButton.setSelectedColor(DnaService.getInstance().getBaseAColor());
    mColorCButton.setSelectedColor(DnaService.getInstance().getBaseCColor());
    mColorGButton.setSelectedColor(DnaService.getInstance().getBaseGColor());
    mColorTButton.setSelectedColor(DnaService.getInstance().getBaseTColor());
    mColorNButton.setSelectedColor(DnaService.getInstance().getBaseNColor());
  }

  @Override
  public final void clicked(ModernClickEvent e) {
    if (e.getMessage().equals(UI.BUTTON_OK)) {
      DnaService.getInstance().setBaseAColor(mColorAButton.getSelectedColor());
      DnaService.getInstance().setBaseCColor(mColorCButton.getSelectedColor());
      DnaService.getInstance().setBaseGColor(mColorGButton.getSelectedColor());
      DnaService.getInstance().setBaseTColor(mColorTButton.getSelectedColor());
      DnaService.getInstance().setBaseNColor(mColorNButton.getSelectedColor());
    }

    super.clicked(e);
  }
}
