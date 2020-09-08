import ROOT
from ROOT import gRandom,TCanvas,TH1F

c1 = TCanvas('c1','Example',200,10,700,500)
hpx = TH1F('hpx','px',100,-4,4)

for i in xrange(25000):
    px = gRandom.Gaus()
    hpx.Fill(px)

hpx.Draw()
c1.Update()

