import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import ij.plugin.*;
import ij.plugin.frame.*;
import java.util.*;
import java.lang.Math;

/*
Adapted from
Olivo-Marin, J.C. 
Extraction of spots in biological images using multiscale products.
Pattern Recognition 35 (2002) 1989-1996.

Bugs on padding pattern were corrected on 03 July 2011
*/

public class Spot_Extraction_Filter implements PlugIn {
    int w, h;
    double kfk0[]={0.0625,0.25,0.375,0.25,0.0625};
    double kfk1[]={0.0625,0,0.25,0,0.375,0,0.25,0,0.0625};
    double kfk2[]={0.0625,0,0,0,0.25,0,0,0,0.375,0,0,0,0.25,0,0,0,0.0625};
    double wts[]={0, 0.55, 0.45};

    public void run(String arg) {
	
	ImagePlus imp=WindowManager.getCurrentImage();

	if (imp==null || imp.getType()==ImagePlus.COLOR_RGB) {
	    IJ.error("This plugin requires a non-RGB image");
	    return;
	}

	int nSlices = imp.getStackSize();

	// Obtaining width and height as a global constants
	w=imp.getWidth();
	h=imp.getHeight();	 

	ImageStack ois = new ImageStack(w, h); // prepare the output imagestack 

	for (int slice=1;slice<=nSlices;slice++) {
	    IJ.showStatus("Processing "+slice+"/"+nSlices+"");
	    IJ.showProgress(slice, nSlices);
	    ImageProcessor ip=imp.getStack().getProcessor(slice).convertToFloat();
	    float[] opix=atrousDecomp((float[])ip.getPixels());
	    ois.addSlice("", opix);
	}
	ImagePlus oimp=new ImagePlus(imp.getShortTitle()+"_ad", ois);
	oimp.show();
	oimp.resetDisplayRange();
    }

    float[] atrousDecomp(float[] pix0) {
	    float[] pix1=convolveKernel(pix0,kfk0);
	    float[] wpix0=subtractImage(pix0,pix1);
	    thresholdSignal(wpix0);
	    powarray(wpix0,wts[0]);

	    float[] pix2=convolveKernel(pix1,kfk1);
	    float[] wpix1=subtractImage(pix1,pix2);
	    thresholdSignal(wpix1);
	    powarray(wpix1,wts[1]);

	    float[] pix3=convolveKernel(pix2,kfk2);
	    float[] wpix2=subtractImage(pix2,pix3);
	    thresholdSignal(wpix2);
	    powarray(wpix2,wts[2]);

	    // Find the cumulative product
	    for (int y=0;y<h;y++) {
		for (int x=0;x<w;x++) {
		    int l=w*y+x;
		    wpix0[l]=wpix1[l]*wpix2[l];
		}
	    }
	    return wpix0;
    }

    float[] convolveKernel(float[] pix, double[] k) {
	int kw=k.length;     // kernel width
	int hkw=(kw-1)/2;    // half kernel width
	int nw=w+kw*2;       // new width after padding
	int nh=h+kw*2;       // new height after padding

	float[] ppix = symmetricalPadding(pix, kw);     // padding to the input image
	double[] gpix=new double[ppix.length];          // temporary array for calculation
	double[] hpix=new double[ppix.length];          // temporary array for calculation
	float[] opix=new float[pix.length];             // array for output

	for (int y=0;y<nh;y++) {
	    for (int x=0;x<nw;x++) {
		int l=nw*y+x;
		if (x<hkw||x>nw-hkw-1)
		    gpix[l]=0;
		else {
		    for (int i=0;i<kw;i++) {
			gpix[l]+=(((double)ppix[nw*y+x+i-hkw])*k[i]);
		    }
		}
	    }
	}

	for (int y=0;y<nh;y++) {
	    for (int x=0;x<nw;x++) {
		int l=nw*y+x;
		if (y<hkw||y>nh-hkw-1)
		    hpix[l]=0;
		else {
		    for (int i=0;i<kw;i++) {
			hpix[l]+=(((double)gpix[nw*(y+i-hkw)+x])*k[i]);
		    }
		}
	    }
	}

	for (int y=0;y<h;y++) {
	    for (int x=0;x<w;x++) {
		opix[w*y+x]=(float)hpix[nw*(kw+y)+kw+x];
	    }
	}
	return opix;
    }

    float[] symmetricalPadding(float[] pix, int p) {
	int nw=w+2*p;        // new width after padding
	int nh=h+2*p;        // new height after padding
	float[] ppix = new float[nw*nh];        // array with padding for output

	for (int y = p; y < h+p; y++) {
	    for (int x = 0; x < p; x++) {
		ppix[nw*y+x]=pix[w*(y-p)+(p-x)];
	    }
	    for (int x = p; x < w+p; x++) {
		ppix[nw*y+x]=pix[w*(y-p)+(x-p)];
	    }
	    for (int x = w+p; x < w+2*p; x++) {
		ppix[nw*y+x]=pix[w*(y-p)+(2*w-x+p-2)];  // modified on 03/07/2011
	    }
	}
	for (int y = 0; y < p; y++) {
	    for (int x = 0; x < p; x++) {
		ppix[nw*y+x]=pix[w*(p-y)+(p-x)];
	    }
	    for (int x = p; x < w+p; x++) {
		ppix[nw*y+x]=pix[w*(p-y)+(x-p)];
	    }
	    for (int x = w+p; x < w+2*p; x++) {
		ppix[nw*y+x]=pix[w*(p-y)+(2*w-x+p-2)];  // modified on 03/07/2011
	    }
	}
	for (int y = h+p; y < h+2*p; y++) {
	    for (int x = 0; x < p; x++) {
		ppix[nw*y+x]=pix[w*(2*h-y+p-2)+(p-x)];  // modified on 03/07/2011
	    }
	    for (int x = p; x < w+p; x++) {
		ppix[nw*y+x]=pix[w*(2*h-y+p-2)+(x-p)];  // modified on 03/07/2011
	    }
	    for (int x = w+p; x < w+2*p; x++) {
		ppix[nw*y+x]=pix[w*(2*h-y+p-2)+(2*w-x+p-2)];  // modified on 03/07/2011
	    }
	}
	return ppix;
    }

    float[] subtractImage(float[] pix1, float[] pix2) {
	int len=pix1.length;
	float[] dpix=new float[len];           // array for output
	for (int i = 0; i < len; i++) {
	    dpix[i]=pix1[i]-pix2[i];
	}
	return dpix;
    }

    void thresholdSignal(float[] pix) {
	float[] ad=new float[pix.length];       // array for temporal store of calculation results.
	double sm=getMedian(pix);               // calculate the median of the input image
	//	IJ.log("sm="+sm+"");

	for (int i=0;i<pix.length;i++) {
	    ad[i]=(float)Math.abs(pix[i]-sm);   // calculate the absolute deviation from the median
	}
	double sig=getMedian(ad);               // Again calculate the median of the above, obtaining "mad"
	//	IJ.log("sig="+sig+"");

	for (int i=0;i<pix.length;i++) {
    	    pix[i]=(pix[i]>=3.0*sig/0.67)? pix[i]:0;       // thresholding
	}
    }

    double getMedian(float[] pix) {
	int len=pix.length;
	float[] a=new float[len];                // temporal array to store the input array for sorting
	System.arraycopy(pix,0,a,0,len);
	Arrays.sort(a);                          // sort the data

	// calculate the median
	double median=0;
	int middle=len/2;
	if (len%2==1)
	    median=a[middle];
	else 
	    median=(a[middle-1]+a[middle])/2.0;
	return median;
    }

    void powarray(float[] a, double b) {
	for (int i=0;i<a.length;i++) {
	    a[i]=(float)Math.pow((double)a[i],b);
	}
    }
}
