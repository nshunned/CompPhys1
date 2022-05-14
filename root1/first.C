// first.C
// This is an unnamed macro.  It essentially functions as a script 
// to specify a series of ROOT commands.
// Usage in ROOT: .X first.C
{
  delete my_hist;
  TH1F *my_hist = new TH1F("my_hist", "My First Histogram;x,y", 100, -5, 5);
  my_hist->FillRandom("gaus", 10000);
  my_hist->Draw();  // the creates a default TCanvas and draws the histogram


  // alternatively, we can make our TCanvas and subdivide it if we want
  auto tc=new TCanvas("tc","My Canvas");  // return pointer to a new TCanvas
  tc->Divide(2,2);   // columns, rows  2x3 TPads in this case
  tc->cd(1);  // TPads are numbered gogin left to right, top to bottom
  my_hist->Draw();
  tc->cd(2);
  my_hist->Draw("E");   // draw with error bars
  tc->cd(3);
  my_hist->Draw("E");   // draw with errorbars
  my_hist->Draw("L,same");   // overlay with lines through the points
  tc->cd(4);
  my_hist->Draw("HBAR");  // draw as horixontal bar chart
  tc->Update()
}
