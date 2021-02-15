
import java.awt.*;
import java.awt.image.*;
import java.io.*;
import javax.swing.*;


public class ImageDisplay {

	JFrame frame;
	JLabel lbIm1;
	JLabel lbIm2;
	BufferedImage imgOne;
	BufferedImage img2;
	int width = 352;
	int height = 288;
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

	/**
	 * Read Image RGB
	 * Reads the image of given width and height at the given imgPath into the provided BufferedImage.
	 */
	private void readImageRGB(int width, int height, String imgPath, BufferedImage img)
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

			int ind = 0;
			for(int y = 0; y < height; y++)
			{
				for(int x = 0; x < width; x++)
				{
					byte a = 0;
					byte r = bytes[ind];
					byte g = bytes[ind+height*width];
					byte b = bytes[ind+height*width*2]; 

					int pix = 0xff000000 | ((r & 0xff) << 16) | ((g & 0xff) << 8) | (b & 0xff);
					// int pix = ((a << 24) + (r << 16) + (g << 8) + b);
					img.setRGB(x,y,pix);
					ind++;
				}
			}
		}
		catch (FileNotFoundException e) 
		{
			e.printStackTrace();
		} 
		catch (IOException e) 
		{
			e.printStackTrace();
		}
	}

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
	private double[][][] convertImgToYUV(int[][][] img, int paramY, int paramU, int paramV)
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
	private int[][][] convertImgToRGB(double[][][] img, int q)
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

	public void showIms(String[] args)
	{
		// Read a parameter from command line
		String imgPath = args[0];
		String paramY = args[1];
		String paramU = args[2];
		String paramV = args[3];
		String paramQ = args[4];
		System.out.println("Image path: " + imgPath);
		System.out.println("The Y parameter was: " + Integer.parseInt(paramY));
		System.out.println("The U parameter was: " + Integer.parseInt(paramU));
		System.out.println("The V parameter was: " + Integer.parseInt(paramV));
		System.out.println("The Q parameter was: " + Integer.parseInt(paramQ));

		// Read in the specified image
		imgOne = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
		readImageRGB(width, height, args[0], imgOne);

		// Processed image
		img2 = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
		byte[] imgBytes = readImg(width, height, imgPath);
		int[][][] rgbImg = convertImgByteToRGB(width, height, imgBytes);
		double[][][] yuvImg = convertImgToYUV(rgbImg, Integer.parseInt(paramY), Integer.parseInt(paramU), Integer.parseInt(paramV));
		int[][][] finalImg = convertImgToRGB(yuvImg, Integer.parseInt(paramQ));
		// calculateColorCount(finalImg);
		readImageIntoBuffer(finalImg, img2);

		// Use label to display the image
		frame = new JFrame();
		GridBagLayout gLayout = new GridBagLayout();
		frame.getContentPane().setLayout(gLayout);

		lbIm1 = new JLabel(new ImageIcon(imgOne));
		lbIm2 = new JLabel(new ImageIcon(img2));

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
