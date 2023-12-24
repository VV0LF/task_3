double psi(double *x, double *par)
{
    return (TMath::Power((2/TMath::Pi()), 0.25) * 1/TMath::Power(par[0], 0.5) * exp(-1 * (TMath::Power(x[0], 2)) / (TMath::Power(par[0], 2))));
}

double U(double x)
{
    if(x < 10 || x > -10) return -0.5;
    return 0;
}

double scalar_expr(double *x, double *par)
{
  double A = par[0];
  TF1 Psi("psi",psi,-10.*A,10.*A,1);
  Psi.SetParameter(0,A);
  // psi* (-ihd/dx)(-ihd/dx) psi 
  // постоянная Планка Math::Hbarcgs()
  double result = - Psi.Eval(x[0])*(Psi.Derivative2(x[0]) + U(x[0]));
  return result;
}

double scal_prod(double *x,double *par)
{
  TF1 scal_exp("sc",scalar_expr,-10.*x[0],10.*x[0],1);
  scal_exp.SetParameter(x[0],1);
  return scal_exp.Integral(-10.*x[0],10.*x[0]);
}

/*
void scalar_prod(double xmin, double xmax)
{
  gStyle->SetLabelSize(0.05,"xy");
  gStyle->SetTitleSize(0.05,"xy");
  TCanvas *s = new TCanvas();
  TF1* scal_exp = new TF1("sp",scalar_expr,xmin,xmax,1);
  scal_exp->Draw();
}
*/

void scalar_prod(double amin, double amax)
{
  gStyle->SetLabelSize(0.05,"xy");
  gStyle->SetTitleSize(0.05,"xy");
  TCanvas *s = new TCanvas();
  TF1 *scalar_pr = new TF1("sp",scal_prod,amin,amax,0);
  scalar_pr->Draw();  
  scalar_pr->SetTitle("<x2> as a function of a");
  scalar_pr->GetXaxis()->SetTitle("a");
  scalar_pr->GetYaxis()->SetTitle("<x2>");
  std::cout << scalar_pr->GetMaximum() << std::endl;
  //s->SaveAs("sc_pr.pdf");
}


void task_3()
{
    Double_t U1 = 0;
    Double_t U2 = -0.5;
    Double_t bord = 10; //A
    scalar_prod(0, 1);
    
}
