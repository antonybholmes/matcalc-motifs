package edu.columbia.rdf.matcalc.toolbox.motifs.app;

import org.jebtk.core.AppVersion;
import org.jebtk.modern.AssetService;
import org.jebtk.modern.help.GuiAppInfo;

public class MotifsInfo extends GuiAppInfo {

  public MotifsInfo() {
    super("Motifs", new AppVersion(3),
        "Copyright (C) 2014-${year} Antony Holmes",
        AssetService.getInstance().loadIcon(MotifsIcon.class, 32),
        AssetService.getInstance().loadIcon(MotifsIcon.class, 128),
        "Search for motifs.");
  }
}
