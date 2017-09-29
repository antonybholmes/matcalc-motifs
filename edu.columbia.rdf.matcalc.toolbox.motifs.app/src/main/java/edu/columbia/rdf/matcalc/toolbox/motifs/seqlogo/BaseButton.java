package edu.columbia.rdf.matcalc.toolbox.motifs.seqlogo;

import java.awt.Color;
import java.awt.Graphics2D;

import org.jebtk.modern.graphics.ImageUtils;
import org.jebtk.modern.graphics.color.ColorPopupMenu;
import org.jebtk.modern.ribbon.RibbonLargeColorSwatchButton;
import org.jebtk.modern.window.ModernWindow;

public class BaseButton extends RibbonLargeColorSwatchButton {
	private static final long serialVersionUID = 1L;

	public static final int SIZE = 24; //Ribbon.COMPACT_ICON_SIZE;

	private String mBase;

	public BaseButton(ModernWindow parent) {
		this(parent, Color.BLACK, "A");
	}

	/**
	 * Instantiates a new color swatch button.
	 *
	 * @param parent the parent
	 * @param color the color
	 */
	public BaseButton(ModernWindow parent, Color color, String base) {
		super(parent, color);

		mBase = base;
	}

	/* (non-Javadoc)
	 * @see org.abh.lib.ui.modern.button.ModernDropDownButton#drawBackground(java.awt.Graphics2D)
	 */
	@Override
	public void drawForegroundAAText(Graphics2D g2) {
		Color color = ((ColorPopupMenu)mMenu).getSelectedColor();

		if (color != null) {
			int x = (getWidth() - 16 - SIZE) / 2; //PADDING;
			int y = (getHeight() - SIZE) / 2;

			Graphics2D g2Temp = ImageUtils.createAAStrokeGraphics(g2);

			try {
				//g2Temp.setColor(mPopup.getSelectedColor());

				//if (isSelected() || mHighlight || mPopupShown) {
				//	g2Temp.setColor(ColorUtils.getTransparentColor25(color));
				//} else {
				//	g2Temp.setColor(color);
				//}

				//g2Temp.fillOval(x, x, SIZE, SIZE);

				/*
				Color c;

				if (isSelected() || mHighlight || mPopupShown) {
					c = ColorUtils.getTransparentColor40(color);
				} else {
					c = color;
				}
				*/

				/*
				getRenderer().fill(g2Temp, 
						color, 
						x,
						y,
						SIZE,
						SIZE);
				*/
				
				g2Temp.setColor(color);
				g2Temp.fillOval(x, y, SIZE, SIZE);

				//g2Temp.setColor(color);
				//g2Temp.drawRect(x, y, SIZE - 1, SIZE - 1);

				/*
				getRenderer().outline(g2Temp, 
						color, 
						x,
						y,
						Ribbon.COMPACT_ICON_SIZE,
						Ribbon.COMPACT_ICON_SIZE);
				*/
			} finally {
				g2Temp.dispose();
			}

			if (color.equals(Color.WHITE)) {
				g2.setColor(Color.BLACK);
			} else {
				g2.setColor(Color.WHITE);
			}

			x += getTextXPosCenter(mBase, SIZE);
			y = getTextYPosCenter(g2, getHeight());

			g2.setFont(BOLD_FONT);
			g2.drawString(mBase, x, y);

			TRIANGLE_ICON.drawIcon(g2, getWidth() - 16, (getHeight() - 16) / 2, 16);
		}
	}
}
