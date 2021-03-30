
import java.awt.*;
import java.awt.image.*;
import java.io.*;
import javax.swing.*;


public class ImageDisplay {

	JFrame frame;
	JLabel lbIm1;
	JLabel lbIm2;
	BufferedImage originalImg;
	BufferedImage processedImg;
	int width = 512;
	int height = 512;
	double [][] RGBtoYUV = {
			{0.299d, 0.587d, 0.114d},
			{0.596d, -0.274d, -0.322d},
			{0.211d, -0.523d, 0.312d}
	};
	double [][] YUVtoRGB = {
			{1.000d, 0.956d, 0.621d},
			{1.000d, -0.272d, -0.647d},
			{1.000d, -1.106d, 1.703d}
	};

	private byte[] readImg(int width, int height, String imgPath)
	{
		try
		{
			int frameLength = width*height*3;

			File file = new File(imgPath);
			RandomAccessFile raf = new RandomAccessFile(file, "r");
			raf.seek(0);

			long len = frameLength;
			byte[] bytes = new byte[(int) len];

			raf.read(bytes);
			return bytes;
		}
		catch (FileNotFoundException e)
		{
			e.printStackTrace();
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
		return new byte[0];
	}

	/**
	 * Read RGB image into buffer
	 * @param img
	 * @param buffer
	 */
	private void readImageIntoBuffer(int[][][] img, BufferedImage buffer)
	{
		int height = img.length; // rows
		int width = img[0].length; // cols
		for(int y = 0; y < height; y++)
		{
			for(int x = 0; x < width; x++)
			{
				byte a = 0;
				byte r = (byte) img[y][x][0];
				byte g = (byte) img[y][x][1];
				byte b = (byte) img[y][x][2];
				int pix = ((a << 24) + (r << 16) + (g << 8) + b);
				buffer.setRGB(x,y,pix);
			}
		}
	}

	/**
	 * convert image bytes to rgb
	 */
	private int[][][] convertImgByteToRGB(int width, int height, byte[] bytes)
	{
		int[][][] img = new int[height][width][3];
		int ind = 0;
		for(int y = 0; y < height; y++)
		{
			for(int x = 0; x < width; x++)
			{
				byte a = 0;
				byte r = bytes[ind];
				byte g = bytes[ind+height*width];
				byte b = bytes[ind+height*width*2];

				img[y][x][0] = (int) r & 0xff;
				img[y][x][1] = (int) g & 0xff;
				img[y][x][2] = (int) b & 0xff;
				ind++;
			}
		}
		return img;
	}

	/**
	 * Convert RGB image to YUV
	 * @param img: RGB image
	 * @return YUV image
	 */
	private double[][][] convertImgRGBToYUV(int[][][] img, int paramY, int paramU, int paramV)
	{
		int height = img.length; // rows
		int width = img[0].length; // cols
		double[][][] yuvImg = new double[height][width][3];
		for(int y = 0; y < height; y++)
		{
			for(int x = 0; x < width; x++)
			{
				int r = img[y][x][0];
				int g = img[y][x][1];
				int b = img[y][x][2];
				double[][] rgb = {
						{(double) r},
						{(double) g},
						{(double) b}
				};
				double[][] yuv = multiplyMatrices(RGBtoYUV, rgb);
				yuvImg[y][x][0] = yuv[0][0];
				yuvImg[y][x][1] = yuv[1][0];
				yuvImg[y][x][2] = yuv[2][0];
			}
		}

		// process lines
		for(int y = 0; y < height; y++)
		{
			double[] yLine = new double[width];
			double[] uLine = new double[width];
			double[] vLine = new double[width];

			// get each channel in a line
			for(int x = 0; x < width; x++)
			{
				yLine[x] = yuvImg[y][x][0];
				uLine[x] = yuvImg[y][x][1];
				vLine[x] = yuvImg[y][x][2];
			}
			// subsampling
			yLine = subsampling(yLine, paramY);
			uLine = subsampling(uLine, paramU);
			vLine = subsampling(vLine, paramV);
			// assign new values
			for(int x = 0; x < width; x++)
			{
				yuvImg[y][x][0] = yLine[x];
				yuvImg[y][x][1] = uLine[x];
				yuvImg[y][x][2] = vLine[x];
			}
		}

		return yuvImg;
	}

	/**
	 * Subsampling a channel with step
	 * @param inputArray array of values in a channel
	 * @param step subsampling step
	 * @return new array with same length
	 */
	private double[] subsampling(double[] inputArray, int step)
	{
		int width = inputArray.length;
		double[] result = new double[width];

		// subsampling
		for(int i = 0; i < width; i++)
		{
			if(i % step == 0)
			{
				result[i] = inputArray[i];
			}
			else
			{
				result[i] = 0;
			}
		}

		// fill the gaps
		for(int i = 0; i < width; i++)
		{
			if(i % step != 0)
			{
				int previous = i - i % step;
				int next = i + step - i % step;
				if(next>=width)
				{
					next = previous;
				}
				result[i] = (result[previous] + result[next]) / 2;
//				System.out.printf("previous %s %s next %s %s%n", previous, result[previous], next, result[next]);
//				result[i] = 0;
			}
		}

		return result;
	}

	/**
	 * Convert YUV image to RGB with Quantization
	 * @param img: YUV image
	 * @param q: quantization
	 * @return RGB image
	 */
	private int[][][] convertImgYUVToRGB(double[][][] img, int q)
	{
		int height = img.length; // rows
		int width = img[0].length; // cols
		int[][][] result = new int[height][width][3];

		for(int y = 0; y < height; y++)
		{
			for(int x = 0; x < width; x++)
			{
				double Y = img[y][x][0];
				double U = img[y][x][1];
				double V = img[y][x][2];
				double[][] yuv = {
						{Y},
						{U},
						{V}
				};
				double[][] rgb = multiplyMatrices(YUVtoRGB, yuv);
//				displayProduct(rgb);
				result[y][x][0] = quantization(rgb[0][0], q);
				result[y][x][1] = quantization(rgb[1][0], q);
				result[y][x][2] = quantization(rgb[2][0], q);
			}
		}

		return result;
	}

	/**
	 * Quantization 8-bit RGB color
	 * @param value
	 * @param q
	 * @return
	 */
	private int quantization(double value, int q)
	{
		// max = 255
		// min = 0
		if(value < 0.0d)
		{
			value = 0.0d;
		}
		if(value > 255.0d)
		{
			value = 255.0d;
		}

		double step = 255.0d / q;
		int newValue = 255;
		double minDiff = 999;
		for (int i=0; i<q; i++)
		{
			int l = (int) (step * i);
			if (Math.abs((double) l- value) < minDiff)
			{
				minDiff = Math.abs((double) l- value);
				newValue = l;
			}
		}

		return newValue;
	}

	/**
	 * convert rgb values to hsv values
	 * @param R
	 * @param G
	 * @param B
	 * @return h:[0, 360] s:[0, 100] v:[0, 100]
	 */
	private double[] convertRGBtoHSV(int R, int G, int B)
	{
		// R, G, B values are divided by 255
		// to change the range from 0..255 to 0..1
		double r = R / 255.0d;
		double g = G / 255.0d;
		double b = B / 255.0d;

		// h, s, v = hue, saturation, value
		double cmax = Math.max(r, Math.max(g, b)); // maximum of r, g, b
		double cmin = Math.min(r, Math.min(g, b)); // minimum of r, g, b
		double diff = cmax - cmin; // diff of cmax and cmin.
		double h = -1, s = -1;
		double v = cmax * 100.0d;

		if (cmax != 0)
		{
			s = (diff / cmax) * 100.0d;
		}
		else
		{
			// s = 0, v is undefined
			s = 0;
			h = -1;
			return new double[]{h, s, v};
		}

		// if cmax and cmax are equal then h = 0
		if (cmax == cmin)
		{
			h = 0;
		}
		else if (cmax == r)
		{
			// between yellow & magenta
			h = (60 * ((g - b) / diff) + 360) % 360;
		}
		else if (cmax == g)
		{
			// between cyan & yellow
			h = (60 * ((b - r) / diff) + 120) % 360;
		}
		else if (cmax == b)
		{
			// between magenta & cyan
			h = (60 * ((r - g) / diff) + 240) % 360;
		}

//		System.out.println("(" + h + " " + s + " " + v + ")");

		return new double[]{h, s, v};
	}

	/**
	 * convert hsv values to rgb values
	 * h:[0, 360] s:[0, 100] v:[0, 100]
	 * @param H
	 * @param S
	 * @param V
	 * @return
	 */
	private int[] convertHSVtoRGB(double H, double S, double V)
	{
		int[] rgb = new int[3];
		double s = S / 100.0d;
		double v = V / 100.0d;
		if(v == 0)
		{
			// achromatic (grey)
			int greyScale = Math.max(0, Math.min((int)Math.round(v * 255), 255));
			rgb[0] = greyScale;
			rgb[1] = greyScale;
			rgb[2] = greyScale;
			return rgb;
		}

		int i;
		double f, p, q, t;
		double r, g, b;
		// sector 0 to 5
		double h = H / 60.0d;
		i = (int)Math.floor(h);
		f = h - i; // factorial part of h
		p = v * ( 1 - s );
		q = v * ( 1 - s * f );
		t = v * ( 1 - s * ( 1 - f ) );
		switch (i) {
			case 0 -> {
				r = v;
				g = t;
				b = p;
			}
			case 1 -> {
				r = q;
				g = v;
				b = p;
			}
			case 2 -> {
				r = p;
				g = v;
				b = t;
			}
			case 3 -> {
				r = p;
				g = q;
				b = v;
			}
			case 4 -> {
				r = t;
				g = p;
				b = v;
			}
			// case 5:
			default -> {
				r = v;
				g = p;
				b = q;
			}
		}
		rgb[0] = Math.max(0, Math.min((int)Math.round(r * 255), 255));
		rgb[1] = Math.max(0, Math.min((int)Math.round(g * 255), 255));
		rgb[2] = Math.max(0, Math.min((int)Math.round(b * 255), 255));

		return rgb;
	}

	private int[] convertRGBtoGreyscale(int R, int G, int B)
	{
		int[] rgb = new int[3];
		int greyScale = Math.max(0, Math.min((int)Math.round(0.299d * R + 0.587d * G + 0.114d * B), 255));
		rgb[0] = greyScale;
		rgb[1] = greyScale;
		rgb[2] = greyScale;
		return rgb;
	}

	private int[][][] filterImg(int[][][] img, double h1, double h2)
	{
		int height = img.length; // rows
		int width = img[0].length; // cols
		int[][][] result = new int[height][width][3];
		for (int y = 0; y < height; y++)
		{
			for (int x = 0; x < width; x++)
			{
				int r = img[y][x][0];
				int g = img[y][x][1];
				int b = img[y][x][2];
				double[] hsv = convertRGBtoHSV(r, g, b);
				double h = hsv[0];
				if (h < h2 && h > h1)
				{
					result[y][x][0] = r;
					result[y][x][1] = g;
					result[y][x][2] = b;
				}
				else
				{
					// System.out.println("(" + h + " " + hsv[1] + " " + hsv[2] + ")");
					// result[y][x] = convertRGBtoGreyscale(r, g, b);
					result[y][x] = convertHSVtoRGB(h, 0.0d, hsv[2]);
				}
			}
		}
		return result;
	}

	public void showIms(String[] args)
	{
		// Read a parameter from command line
		String imgPath = args[0];
		String paramH1 = args[1];
		String paramH2 = args[2];
		System.out.println("Image path: " + imgPath);
		System.out.println("The h1 parameter was: " + Integer.parseInt(paramH1));
		System.out.println("The h2 parameter was: " + Integer.parseInt(paramH2));

		// Read in the original image
		originalImg = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
		byte[] originalImgBytes = readImg(width, height, imgPath);
		int[][][] originalRGBImg = convertImgByteToRGB(width, height, originalImgBytes);
		readImageIntoBuffer(originalRGBImg, originalImg);

		// Processed image
		processedImg = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
		byte[] imgBytes = readImg(width, height, imgPath);
		int[][][] rgbImg = convertImgByteToRGB(width, height, imgBytes);
		int[][][] finalImg = filterImg(rgbImg, Double.parseDouble(paramH1), Double.parseDouble(paramH2));
		readImageIntoBuffer(finalImg, processedImg);

		// Use label to display the image
		frame = new JFrame();
		GridBagLayout gLayout = new GridBagLayout();
		frame.getContentPane().setLayout(gLayout);

		lbIm1 = new JLabel(new ImageIcon(originalImg));
		lbIm2 = new JLabel(new ImageIcon(processedImg));

		GridBagConstraints c = new GridBagConstraints();
		c.fill = GridBagConstraints.HORIZONTAL;
		c.anchor = GridBagConstraints.CENTER;
		c.weightx = 0.5;
		c.gridx = 0;
		c.gridy = 0;

		c.fill = GridBagConstraints.HORIZONTAL;
		c.gridx = 0;
		c.gridy = 1;
		frame.getContentPane().add(lbIm1);
		frame.getContentPane().add(lbIm2);

		frame.pack();
		frame.setVisible(true);
	}

	public static double[][] multiplyMatrices(double[][] firstMatrix, double[][] secondMatrix) {
		double[][] result = new double[firstMatrix.length][secondMatrix[0].length];

		for (int row = 0; row < result.length; row++) {
			for (int col = 0; col < result[row].length; col++) {
				result[row][col] = multiplyMatricesCell(firstMatrix, secondMatrix, row, col);
			}
		}

		return result;
	}

	public static double multiplyMatricesCell(double[][] firstMatrix, double[][] secondMatrix, int row, int col) {
		double cell = 0;
		for (int i = 0; i < secondMatrix.length; i++) {
			cell += firstMatrix[row][i] * secondMatrix[i][col];
		}
		return cell;
	}

	public static void displayProduct(double[][] product) {
		System.out.println("Product of two matrices is: ");
		for(double[] row : product) {
			for (double column : row) {
				System.out.print(column + "    ");
			}
			System.out.println();
		}
	}

	public static void main(String[] args) {
		ImageDisplay ren = new ImageDisplay();
		ren.showIms(args);
	}

}
