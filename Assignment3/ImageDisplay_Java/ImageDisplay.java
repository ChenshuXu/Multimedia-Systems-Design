
import java.awt.*;
import java.awt.image.*;
import java.io.*;
import java.util.Arrays;
import javax.swing.*;


public class ImageDisplay {

	JFrame frame;
	JLabel lbIm1;
	JLabel lbIm2;
	JLabel lbIm3;
	BufferedImage originalImg;
	BufferedImage processedImg;
	BufferedImage processedImg2;
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


	private void printDoubleMatrix(double[][] matrix)
	{
		int row = matrix.length;
		int col = matrix[0].length;
		int r, c;
		for (r = 0; r < row; r++)
		{
			for (c = 0; c < col; c++) {
				System.out.printf("%f\t", matrix[r][c]);
			}
			System.out.println();
		}
	}

	private void printIntMatrix(int[][] matrix)
	{
		int row = matrix.length;
		int col = matrix[0].length;
		int r, c;
		for (r = 0; r < row; r++)
		{
			for (c = 0; c < col; c++) {
				System.out.printf("%d\t", matrix[r][c]);
			}
			System.out.println();
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


	/**
	 * convert a rgb pixel to grey scale pixel
	 * @param R
	 * @param G
	 * @param B
	 * @return
	 */
	private int[] convertRGBtoGreyscale(int R, int G, int B)
	{
		int[] rgb = new int[3];
		int greyScale = Math.max(0, Math.min((int)Math.round(0.299d * R + 0.587d * G + 0.114d * B), 255));
		rgb[0] = greyScale;
		rgb[1] = greyScale;
		rgb[2] = greyScale;
		return rgb;
	}

	/**
	 * filter out pixels with hue value between h1 and h2
	 * @param img input image
	 * @param h1 hue value 1
	 * @param h2 hue value 2
	 * @return filter image
	 */
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

	/**
	 * 8 * 8 dct
	 * @param imageValues a 8 by 8 matrix
	 * @return a 8 by 8 matrix result
	 */
	private double[][] dctTransform(int[][] imageValues)
	{
		int u, v, x, y;
		double[][] dct = new double[8][8];
		double cu, cv, dct1, sum;

		for (u = 0; u < 8; u++)
		{
			for (v = 0; v < 8; v++)
			{
				// get cu amd cv
				if (u == 0)
				{
					cu = 1 / Math.sqrt(2);
				}
				else
				{
					cu = 1;
				}

				if (v == 0)
				{
					cv = 1 / Math.sqrt(2);
				}
				else
				{
					cv = 1;
				}

				// sum of cosine
				sum = 0;
				for (x = 0; x < 8; x++)
				{
					for (y = 0; y < 8; y++)
					{
						dct1 = imageValues[x][y] * Math.cos((2 * x + 1.0d) * u * Math.PI / 16.0d) *
								Math.cos((2 * y + 1.0d) * v * Math.PI / 16.0d);
						sum += dct1;
					}
				}
				dct[u][v] = cu * cv * sum / 4.0d;
			}
		}

		return dct;
	}


	/**
	 * inverse dct on 8 by 8 block
	 * @param dctValues 8 by 8 block
	 * @return same size
	 */
	private double[][] inverseDctTransform(double[][] dctValues)
	{
		int u, v, x, y;
		double[][] imageValues = new double[8][8];
		double cu, cv, idct1, sum;

		for (x = 0; x < 8; x++)
		{
			for (y = 0; y < 8; y++)
			{
				// sum of cosine
				sum = 0;
				for (u = 0; u< 8; u++)
				{
					for (v = 0; v < 8; v++)
					{
						// get cu amd cv
						if (u == 0)
						{
							cu = 1 / Math.sqrt(2);
						}
						else
						{
							cu = 1;
						}

						if (v == 0)
						{
							cv = 1 / Math.sqrt(2);
						}
						else
						{
							cv = 1;
						}
						idct1 = cu * cv * dctValues[u][v] *
								Math.cos((2 * x + 1.0d) * u * Math.PI / 16.0d) *
								Math.cos((2 * y + 1.0d) * v * Math.PI / 16.0d);
						sum += idct1;
					}
				}
				imageValues[x][y] = sum / 4.0d;
			}
		}

		return imageValues;
	}


	/**
	 * get the first m coefficients in zig-zag order, zeroing out others
	 * @param matrix a 8 by 8 matrix
	 * @param limit first m coefficients to keep
	 * @return
	 */
	private double[][] getFirstMCoefficients(double[][] matrix, int limit)
	{
		boolean overLimit = false;

		// Variables to track the size of the matrix
		int N = 8;
		int M = 8;

		int row = 0, column = 0;

		int direction = 1;

		double[][] result = new double[8][8];
		int count = 0;

		while (row < 8 && column < 8)
		{
			// check limit
			if (overLimit)
			{
				result[row][column] = 0.0d;
			}
			else {
				result[row][column] = matrix[row][column];
				if (count >= limit-1)
				{
					overLimit = true;
				}
			}
			count += 1;

			// move along in the current diagonal
			// if going up [i, j] -> [i - 1, j + 1]
			// if going down [i, j] -> [i + 1][j - 1]
			int new_row = row + (direction == 1 ? -1 : 1);
			int new_column = column + (direction == 1 ? 1 : -1);

			// Checking if the next element in the diagonal is within the
			// bounds of the matrix or not. If it's not within the bounds,
			// we have to find the next head.
			if (new_row < 0 || new_row == 8 || new_column < 0 || new_column == 8) {

				// If the current diagonal was going in the upwards
				// direction.
				if (direction == 1) {

					// For an upwards going diagonal having [i, j] as its tail
					// If [i, j + 1] is within bounds, then it becomes
					// the next head. Otherwise, the element directly below
					// i.e. the element [i + 1, j] becomes the next head
					if (column == 8-1)
					{
						row += 1;
					}
					if (column < 8-1)
					{
						column += 1;
					}

				} else {

					// For a downwards going diagonal having [i, j] as its tail
					// if [i + 1, j] is within bounds, then it becomes
					// the next head. Otherwise, the element directly below
					// i.e. the element [i, j + 1] becomes the next head
					if (row == 8-1)
					{
						column += 1;
					}
					if (row < 8-1)
					{
						row += 1;
					}
				}

				// flip the direction
				direction = 1 - direction;

			} else {
				row = new_row;
				column = new_column;
			}
		}
		return result;
	}


	/**
	 * use dct compression to compress a 8 by 8 block
	 * @param block image block
	 * @param limit  limit of coefficients
	 * @return compressed image
	 */
	private int[][] dctCompressionBlock(int[][] block, int limit)
	{
		double[][] dctValues = dctTransform(block);
		double[][] partDctValues = getFirstMCoefficients(dctValues, limit);
		double[][] idctValues = inverseDctTransform(partDctValues);

		// round values
		int[][] imageValues = new int[8][8];
		int r, c;
		for (r = 0; r < 8; r++)
		{
			for (c = 0; c < 8; c++)
			{
				int v = (int)Math.round(idctValues[r][c]);
				if (v > 255){
					v = 255;
				}
				if (v < 0){
					v = 0;
				}
				imageValues[r][c] = v;
			}
		}
		return imageValues;
	}


	/**
	 * use dct compression to compress a channel
	 * @param channel
	 * @param numberOfCoefficients
	 * @return
	 */
	private int[][] dctCompressionChannel(int[][] channel, int numberOfCoefficients)
	{
		int m = (int) Math.ceil(((double)numberOfCoefficients)/4096.0d);

		// 8 by 8 block
		int row = channel.length;
		int col = channel[0].length;
		int[][] result = new int[row][col];
		int r, c;
		int[][] block = new int[8][8];
		for (r = 0; r < row; r += 8)
		{
			for (c = 0; c < col; c += 8)
			{
				int i, j; // index in block
				// get block
				for (i = 0; i < 8; i++)
				{
					for (j = 0; j < 8; j++)
					{
						block[i][j] = channel[r + i][c + j];
					}
				}
				// process block
				block = dctCompressionBlock(block, m);
//				printIntMatrix(block);
				// put block on channel
				for (i = 0; i < 8; i++)
				{
					for (j = 0; j < 8; j++)
					{
						result[r + i][c + j] = block[i][j];
					}
				}
			}
		}
		return result;
	}


	/**
	 * use dct compression to compress image
	 * @param rgbImage
	 * @param numberOfCoefficients specific: 262144, 65536, 16384,4096, 1024, 256, 64, 16, 4, 1
	 * @return
	 */
	private int[][][] dctCompressionImage(int[][][] rgbImage, int numberOfCoefficients)
	{
		int row = rgbImage.length;
		int col = rgbImage[0].length;
		int[][][] result = new int[row][col][3];
		// split image to three channels
		int[][] red = new int[row][col];
		int[][] green = new int[row][col];
		int[][] blue = new int[row][col];
		int r, c;
		for (r = 0; r < row; r++)
		{
			for (c = 0; c < col; c++)
			{
				red[r][c] = rgbImage[r][c][0];
				green[r][c] = rgbImage[r][c][1];
				blue[r][c] = rgbImage[r][c][2];
			}
		}

		red = dctCompressionChannel(red, numberOfCoefficients);
		green = dctCompressionChannel(green, numberOfCoefficients);
		blue = dctCompressionChannel(blue, numberOfCoefficients);

		// merge three channels to image
		for (r = 0; r < row; r++)
		{
			for (c = 0; c < col; c++)
			{
				result[r][c][0] = red[r][c];
				result[r][c][1] = green[r][c];
				result[r][c][2] = blue[r][c];
			}
		}
		return result;
	}


	/**
	 * apply low pass and high pass on an array, dwt
	 * @param array length must be divided by 2
	 * @return same size as input array
	 */
	private double[] dwtOnArray(double[] array)
	{
		int len = array.length;
		int halfLen = len / 2;
		double[] result = new double[len];
		int i;
		for (i = 0; i < len; i += 2)
		{
			result[i/2] = (array[i] + array[i+1]) / 2.0d;
			result[halfLen + i/2] = (array[i] - array[i+1]) / 2.0d;
		}
		return result;
	}


	/**
	 * apply inverse dwt on dwt array
	 * @param dwtArray length must be divided by 2
	 * @return same size as input array
	 */
	private double[] inverseDwtOnArray(double[] dwtArray)
	{
		int len = dwtArray.length;
		int halfLen = len / 2;
		double[] result = new double[len];
		int i;
		for (i = 0; i < len; i += 2)
		{
			double sum = dwtArray[i/2]*2.0d;
			double diff = dwtArray[halfLen + i/2]*2.0d;
			result[i] = (sum + diff) / 2.0d;
			result[i+1] = (sum - diff) / 2.0d;
		}
		return result;
	}


	/**
	 * apply dwt compression on channel pixels
	 * @param channel pixels in a color channel
	 * @param numberOfCoefficients can be divided by 4,
	 *                   specific: 262144, 65536, 16384,4096, 1024, 256, 64, 16, 4, 1
	 * @return channel pixels after compression
	 */
	private int[][] dwtCompressionChannel(int[][] channel, int numberOfCoefficients)
	{
		int row = channel.length;
		int col = channel[0].length;
		int r, c;
		int[][] result = new int[row][col];

		int totalLevel = (int)(Math.log(row) / Math.log(2)); // 9 for 512 * 512 image

		double[][] dwtChannelValues = new double[row][col];
		// convert to double channel values
		for (r=0; r<row; r++)
		{
			for (c=0; c<col; c++)
			{
				dwtChannelValues[r][c] = channel[r][c];
			}
		}

		int gridRow, gridCol;
		double[] currentRow, currentCol;
		double[] newRow, newCol;
		// apply dwt on channel
		for (int l = 0; l < totalLevel; l++)
		{
			gridCol = col / (int)(Math.pow(2, l));
			gridRow = row / (int)(Math.pow(2, l));
			// apply inverse dwt on columns
			for (c = 0; c < gridCol; c++)
			{
				// get current column
				currentCol = new double[gridRow];
				for (r = 0; r < gridRow; r++)
				{
					currentCol[r] = dwtChannelValues[r][c];
				}
				newCol = dwtOnArray(currentCol);
				// update channel values
				for (r = 0; r < gridRow; r++)
				{
					dwtChannelValues[r][c] = newCol[r];
				}
			}
			// apply inverse dwt on rows
			for (r = 0; r < gridRow; r++)
			{
				// get current row
				currentRow = new double[gridCol];
				for (c = 0; c < gridCol; c++)
				{
					currentRow[c] = dwtChannelValues[r][c];
				}
				newRow = dwtOnArray(currentRow);
				// update channel values
				for (c = 0; c < gridCol; c++)
				{
					dwtChannelValues[r][c] = newRow[c];
				}
			}
		}

		// keep the values in range
		int width = (int)Math.round(Math.sqrt(numberOfCoefficients));
		for (r = 0; r < row; r++)
		{
			for (c = 0; c < col; c++)
			{
				if (r > width || c > width)
				{
					dwtChannelValues[r][c] = 0.0d;
				}
			}
		}

		// inverse dwt
		for (int l = totalLevel-1; l >= 0; l--)
		{
			gridCol = col / (int)(Math.pow(2, l));
			gridRow = row / (int)(Math.pow(2, l));
			// apply inverse dwt on columns
			for (c = 0; c < gridCol; c++)
			{
				// get current column
				currentCol = new double[gridRow];
				for (r = 0; r < gridRow; r++)
				{
					currentCol[r] = dwtChannelValues[r][c];
				}
				newCol = inverseDwtOnArray(currentCol);
				// update channel values
				for (r = 0; r < gridRow; r++)
				{
					dwtChannelValues[r][c] = newCol[r];
				}
			}
			// apply inverse dwt on rows
			for (r = 0; r < gridRow; r++)
			{
				// get current row
				currentRow = new double[gridCol];
				for (c = 0; c < gridCol; c++)
				{
					currentRow[c] = dwtChannelValues[r][c];
				}
				newRow = inverseDwtOnArray(currentRow);
				// update channel values
				for (c = 0; c < gridCol; c++)
				{
					dwtChannelValues[r][c] = newRow[c];
				}
			}
		}

		// convert to int
		for (r = 0; r < row; r++)
		{
			for (c = 0; c < col; c++)
			{
				int v = (int)Math.round(dwtChannelValues[r][c]);
				if (v > 255){
					v = 255;
				}
				if (v < 0)
				{
					v = 0;
				}
				result[r][c] = v;
			}
		}

		return result;
	}


	private int[][][] dwtCompressionImage(int[][][] rgbImage, int numberOfCoefficients)
	{
		int row = rgbImage.length;
		int col = rgbImage[0].length;
		int[][][] result = new int[row][col][3];
		// split image to three channels
		int[][] red = new int[row][col];
		int[][] green = new int[row][col];
		int[][] blue = new int[row][col];
		int r, c;
		for (r = 0; r < row; r++)
		{
			for (c = 0; c < col; c++)
			{
				red[r][c] = rgbImage[r][c][0];
				green[r][c] = rgbImage[r][c][1];
				blue[r][c] = rgbImage[r][c][2];
			}
		}

		red = dwtCompressionChannel(red, numberOfCoefficients);
		green = dwtCompressionChannel(green, numberOfCoefficients);
		blue = dwtCompressionChannel(blue, numberOfCoefficients);

		// merge three channels to image
		for (r = 0; r < row; r++)
		{
			for (c = 0; c < col; c++)
			{
				result[r][c][0] = red[r][c];
				result[r][c][1] = green[r][c];
				result[r][c][2] = blue[r][c];
			}
		}
		return result;
	}


	public void showIms(String[] args)
	{
		// Read a parameter from command line
		String imgPath = args[0];
		String param1 = args[1];
		System.out.println("Image path: " + imgPath);
		System.out.println("The h1 parameter was: " + Integer.parseInt(param1));

		// Read in the original image
		originalImg = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
		byte[] originalImgBytes = readImg(width, height, imgPath);
		int[][][] originalRGBImg = convertImgByteToRGB(width, height, originalImgBytes);
		readImageIntoBuffer(originalRGBImg, originalImg);

		// Processed image
		processedImg = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
		processedImg2 = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
		byte[] imgBytes = readImg(width, height, imgPath);
		int[][][] rgbImg = convertImgByteToRGB(width, height, imgBytes);

		// test code here
//		double[] test = new double[] {9., 7., 3., 5., 3., 5., 7., 9., 10., 12.,};
//		double[] lowHighPass = dwtOnArray(test);
//		for (int i = 0; i < test.length; i++)
//		{
//			System.out.printf("%f\t", lowHighPass[i]);
//		}
//		System.out.print("\n");
//
//		double[] inverse = inverseDwtOnArray(lowHighPass);
//		for (int i = 0; i < test.length; i++)
//		{
//			System.out.printf("%f\t", inverse[i]);
//		}

		// add code here
		int numberOfCoefficients = Integer.parseInt(param1);
		int[][][] dctCompressedImage = dctCompressionImage(rgbImg, numberOfCoefficients);
		readImageIntoBuffer(dctCompressedImage, processedImg);
		int[][][] dwtCompressedImage = dwtCompressionImage(rgbImg, numberOfCoefficients);
		readImageIntoBuffer(dwtCompressedImage, processedImg2);
		// end of code here

		// Use label to display the image
		frame = new JFrame();
		GridBagLayout gLayout = new GridBagLayout();
		frame.getContentPane().setLayout(gLayout);

		lbIm1 = new JLabel(new ImageIcon(originalImg));
		lbIm2 = new JLabel(new ImageIcon(processedImg));
		lbIm3 = new JLabel(new ImageIcon(processedImg2));

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
		frame.getContentPane().add(lbIm3);

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
