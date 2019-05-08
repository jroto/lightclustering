#ifndef __MYSTYLE_C
#define __MYSTYLE_C


#include <iostream>

#include "TCanvas.h"
//#include "TRoot.h"
#include "TStyle.h"

using namespace std;

TStyle* MyStyle(bool b_2D=false);

void SetMyStyle (bool b_2D=false)
{
  std::cout << "\nApplying My style settings...\n" << std::endl ;
  //TStyle* MyStyle = 
  MyStyle(b_2D);
  gROOT->SetStyle("MyStyle");
  gROOT->ForceStyle();
}

TStyle* MyStyle(bool b_2D) 
{
  TStyle *MyStyle = new TStyle("MyStyle","My style");

  // use plain black on white colors
  Int_t icol=0; // WHITE
  MyStyle->SetFrameBorderMode(icol);
  MyStyle->SetFrameFillColor(icol);
  MyStyle->SetCanvasBorderMode(icol);
  MyStyle->SetCanvasColor(icol);
  MyStyle->SetPadBorderMode(icol);
  MyStyle->SetPadColor(icol);
  MyStyle->SetStatColor(icol);
  //MyStyle->SetFillColor(icol); // don't use: white fill color floa *all* objects

  // set the paper & margin sizes
  MyStyle->SetPaperSize(20,26);

  // set margin sizes
  MyStyle->SetPadTopMargin(0.05);
  if (b_2D) MyStyle->SetPadRightMargin(0.16);
  else MyStyle->SetPadRightMargin(0.05);
  MyStyle->SetPadBottomMargin(0.16);
  MyStyle->SetPadLeftMargin(0.16);

  // set title offsets (for axis label)
  MyStyle->SetTitleXOffset(1.25);
  MyStyle->SetTitleYOffset(1.25); // was 1.25

  // use large fonts
  //Int_t font=72; // Helvetica italics
  Int_t font=42; // Helvetica
  Double_t tsize=0.06;
  MyStyle->SetTextFont(font);

  MyStyle->SetTextSize(tsize);
  MyStyle->SetLabelFont(font,"x");
  MyStyle->SetTitleFont(font,"x");
  MyStyle->SetLabelFont(font,"y");
  MyStyle->SetTitleFont(font,"y");
  MyStyle->SetLabelFont(font,"z");
  MyStyle->SetTitleFont(font,"z");
  
  MyStyle->SetLabelSize(tsize,"x");
  MyStyle->SetTitleSize(tsize,"x");
  MyStyle->SetLabelSize(tsize,"y");
  MyStyle->SetTitleSize(tsize,"y");
  MyStyle->SetLabelSize(tsize,"z");
  MyStyle->SetTitleSize(tsize,"z");

MyStyle->SetTitleAlign();
MyStyle->SetTitleSize(0.075,"a");
MyStyle->SetTitleX(0.5);
MyStyle->SetTitleY(0.9);
  // use bold lines and markers
//  MyStyle->SetMarkerStyle(20);
//  MyStyle->SetMarkerSize(1.2);
//  MyStyle->SetHistLineWidth(2);
  MyStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  // get rid of X error bars and y error bar caps
  //MyStyle->SetErrorX(0.001);

  // do not display any of the standard histogram decorations
  MyStyle->SetOptTitle(1);
  //MyStyle->SetOptStat(1111);
  MyStyle->SetOptStat(0);
  //MyStyle->SetOptFit(1111);
  MyStyle->SetOptFit(0);

  // put tick marks on top and RHS of plots
  MyStyle->SetPadTickX(1);
  MyStyle->SetPadTickY(1);

  return MyStyle;

}

#endif // __MYSTYLE_C
