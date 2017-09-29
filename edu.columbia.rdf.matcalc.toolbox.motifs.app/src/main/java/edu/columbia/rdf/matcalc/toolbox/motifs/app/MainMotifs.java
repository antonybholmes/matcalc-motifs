package edu.columbia.rdf.matcalc.toolbox.motifs.app;

import java.awt.FontFormatException;
import java.io.IOException;

import javax.swing.UnsupportedLookAndFeelException;
import javax.xml.parsers.ParserConfigurationException;

import org.jebtk.core.AppService;
import org.jebtk.modern.ColorTheme;
import org.jebtk.modern.help.GuiAppInfo;
import org.jebtk.modern.theme.ThemeService;
import edu.columbia.rdf.matcalc.MainMatCalc;
import edu.columbia.rdf.matcalc.ModuleLoader;

import org.xml.sax.SAXException;

import edu.columbia.rdf.matcalc.bio.BioModuleLoader;
import edu.columbia.rdf.matcalc.toolbox.motifs.MotifsModule;




public class MainMotifs {
	public static final void main(String[] args) throws FontFormatException, IOException, SAXException, ParserConfigurationException, ClassNotFoundException, InstantiationException, IllegalAccessException, UnsupportedLookAndFeelException {
		AppService.getInstance().setAppInfo("motifs");
		
		ThemeService.getInstance().setTheme(ColorTheme.PURPLE);
		
		
		GuiAppInfo info = new MotifsInfo();
		
		
		ModuleLoader ml = new BioModuleLoader().addModule(MotifsModule.class);
		
		MainMatCalc.main(info, ml);
	}
}
