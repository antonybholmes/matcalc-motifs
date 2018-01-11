package edu.columbia.rdf.matcalc.toolbox.motifs.seqlogo;

import org.jebtk.core.AppVersion;
import org.jebtk.modern.UIService;
import org.jebtk.modern.help.GuiAppInfo;

public class SeqLogoInfo extends GuiAppInfo {

  public SeqLogoInfo() {
    super("SeqLogo", new AppVersion(5), "Copyright (C) 2016 Antony Holmes",
        UIService.getInstance().loadIcon(SeqLogoIcon.class, 32),
        UIService.getInstance().loadIcon(SeqLogoIcon.class, 128),
        "View motifs graphically.");
  }

}
