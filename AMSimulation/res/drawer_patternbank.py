#!/usr/bin/env python

from rootdrawing import *
from roothelper import *
from array import array

# ______________________________________________________________________________
# Configurations

sections = {}
sections["fixed"     ] = True
sections["projective"] = False

drawerInit = DrawerInit()
gStyle.SetPadRightMargin(0.1)

EOS = "/eos/uscms/store/user/l1upgrades/SLHC/GEN/620_SLHC12p1_results/"

imgdir = "figures_patternbank/"

chain  = TChain("patternBank", "")
chain2 = TChain("patternBankStats", "")

if not imgdir.endswith("/"):  imgdir += "/"
if gSystem.AccessPathName(imgdir):
    gSystem.mkdir(imgdir)


# ______________________________________________________________________________
# Fixed

if sections["fixed"]:

    gStyle.SetTitleSize(0.05, "Y")
    latex.SetTextSize(0.05)

    in_superstrips = ["ss32", "ss64", "ss128", "ss256", "ss512", "ss1024"]
    #in_superstrips = ["ss256"]

    def bookCoverage(bank):
        superstrips = []
        graphs = {}

        for ss in in_superstrips:
            chain.Reset(); chain2.Reset()
            infile = (EOS + "/" + bank) % ss

            # Pattern bank statistics
            chain2.Add(infile)
            assert(chain2.GetEntries() > 0)
            chain2.GetEntry(0)
            coverage = chain2.coverage
            count = chain2.count

            # Pattern bank
            chain.Add(infile)
            assert(chain.GetEntries() > 0)
            chain.SetBranchStatus("superstripIds", 0)
            npatterns = chain.GetEntriesFast()
            npoints = 200
            every = npatterns / npoints

            xvalues, yvalues = [], []
            x, x_0p9, integral = 0, 0, 0
            for i in xrange(npatterns):
                chain.GetEntry(i)
                frequency = chain.frequency
                if i == x:
                    y = float(integral) / float(count) * coverage
                    xvalues.append(i); yvalues.append(y)
                    print "..", i, y
                    if y < 0.9:  x_0p9 = x
                    x += every

                elif i == npatterns - 1:
                    y = float(integral) / float(count) * coverage
                    xvalues.append(i); yvalues.append(y)
                    print "..", i, y
                    if ((integral + frequency) != count):
                        print "ERROR: ", integral, frequency, count
                    #assert((integral + frequency) == count)

                integral += frequency

            superstrips.append((ss, npatterns, coverage, x_0p9))
            print superstrips[-1]
            npoints = len(xvalues)
            gname = "gr_%s" % ss
            gr = TGraph(npoints, array('d', xvalues), array('d', yvalues))
            gr.SetName(gname)
            graphs[gname] = gr
        return (superstrips, graphs)

    def drawCoverage(superstrips, graphs, xmin=0, xmax=1e8, tower="tt27"):
        hframe = TH1F("hframe", "; # of patterns; running estimate for coverage", 100, xmin, xmax)
        hframe.SetStats(0); hframe.SetMinimum(0); hframe.SetMaximum(1.2)
        hframe.SetNdivisions(510, "Y")

        # Style
        for i, ss in enumerate(superstrips):
            gr = graphs["gr_%s" % ss[0]]
            gr.SetLineWidth(2); gr.SetLineStyle(1); gr.SetMarkerSize(0)
            gr.SetLineColor(paletteSet1[i])

        # Legend
        moveLegend(0.66,0.15,0.96,0.45); legend.Clear()
        for i, ss in enumerate(superstrips):
            gr = graphs["gr_%s" % ss[0]]
            legend.AddEntry(gr, ss[0], "l")

        # Draw
        hframe.Draw()
        for y in [0.5, 0.8, 0.9, 0.95, 1.0]:
            line.DrawLine(xmin, y, xmax, y)
        for i, ss in enumerate(superstrips):
            gr = graphs["gr_%s" % ss[0]]
            gr.Draw("C")
        legend.Draw()
        CMS_label()
        save(imgdir, "fixed_sorted_%s" % tower, dot_root=True)

        # Zoom in
        hframe.Draw()
        hframe.GetXaxis().SetRangeUser(0, xmax/50)
        for y in [0.5, 0.8, 0.9, 0.95, 1.0]:
            line.DrawLine(xmin, y, xmax/50, y)
        for i, ss in enumerate(superstrips):
            gr = graphs["gr_%s" % ss[0]]
            gr.Draw("C")
        legend.Draw()
        CMS_label()
        save(imgdir, "fixed_sorted_zoom_%s" % tower)

        # Print out
        for i, ss in enumerate(superstrips):
            print '{0:7} {1:10d}  {2:5.4f}'.format(ss[0], ss[1], ss[2])
        print
        for i, ss in enumerate(superstrips):
            print '{0:7} {1:10d}  {2:5.4f}'.format(ss[0], ss[3], 0.9)

        donotdelete = [hframe]
        return donotdelete

    # Barrel 2 GeV
    bank = "patternBank_sp16_%s_tt27_400M.root"
    (superstrips, graphs) = bookCoverage(bank)
    d = drawCoverage(superstrips, graphs, xmax=5e7, tower="tt27")

    # Barrel 3 GeV
    #bank = "patternBank_sp16_%s_tt27_pt3_400M.root"
    #(superstrips, graphs) = bookCoverage(bank)
    #d = drawCoverage(superstrips, graphs, xmax=5e7, tower="tt27_pt3")


# ______________________________________________________________________________
# Projective

if sections["projective"]:

    gStyle.SetTitleSize(0.05, "Y")
    latex.SetTextSize(0.05)

    in_superstrips = ["600x0", "400x0", "200x0", "200x1", "100x2", "20x10"]
    #in_superstrips = ["400x0"]

    def bookCoverage(bank):
        superstrips = []
        graphs = {}

        for ss in in_superstrips:
            chain.Reset(); chain2.Reset()
            infile = (EOS + "/" + bank) % ss

            # Pattern bank statistics
            chain2.Add(infile)
            assert(chain2.GetEntries() > 0)
            chain2.GetEntry(0)
            coverage = chain2.coverage
            count = chain2.count

            # Pattern bank
            chain.Add(infile)
            assert(chain.GetEntries() > 0)
            chain.SetBranchStatus("superstripIds", 0)
            npatterns = chain.GetEntriesFast()
            npoints = 200
            every = npatterns / npoints

            xvalues, yvalues = [], []
            x, x_0p9, integral = 0, 0, 0
            for i in xrange(npatterns):
                chain.GetEntry(i)
                frequency = chain.frequency
                if i == x:
                    y = float(integral) / float(count) * coverage
                    xvalues.append(i); yvalues.append(y)
                    print "..", i, y
                    if y < 0.9:  x_0p9 = x
                    x += every

                elif i == npatterns - 1:
                    y = float(integral) / float(count) * coverage
                    xvalues.append(i); yvalues.append(y)
                    print "..", i, y
                    if ((integral + frequency) != count):
                        print "ERROR: ", integral, frequency, count
                    #assert((integral + frequency) == count)

                integral += frequency

            superstrips.append((ss, npatterns, coverage, x_0p9))
            print superstrips[-1]
            npoints = len(xvalues)
            gname = "gr_%s" % ss
            gr = TGraph(npoints, array('d', xvalues), array('d', yvalues))
            gr.SetName(gname)
            graphs[gname] = gr
        return (superstrips, graphs)

    def drawCoverage(superstrips, graphs, xmin=0, xmax=1e8, tower="tt27"):
        hframe = TH1F("hframe", "; # of patterns; running estimate for coverage", 100, xmin, xmax)
        hframe.SetStats(0); hframe.SetMinimum(0); hframe.SetMaximum(1.2)
        hframe.SetNdivisions(510, "Y")

        # Style
        for i, ss in enumerate(superstrips):
            gr = graphs["gr_%s" % ss[0]]
            gr.SetLineWidth(2); gr.SetLineStyle(1); gr.SetMarkerSize(0)
            gr.SetLineColor(paletteSet1[i])

        # Legend
        moveLegend(0.66,0.15,0.96,0.45); legend.Clear()
        for i, ss in enumerate(superstrips):
            gr = graphs["gr_%s" % ss[0]]
            legend.AddEntry(gr, ss[0], "l")

        # Draw
        hframe.Draw()
        for y in [0.5, 0.8, 0.9, 0.95, 1.0]:
            line.DrawLine(xmin, y, xmax, y)
        for i, ss in enumerate(superstrips):
            gr = graphs["gr_%s" % ss[0]]
            gr.Draw("C")
        legend.Draw()
        CMS_label()
        save(imgdir, "projective_sorted_%s" % tower, dot_root=True)

        # Zoom in
        hframe.Draw()
        hframe.GetXaxis().SetRangeUser(0, xmax/50)
        for y in [0.5, 0.8, 0.9, 0.95, 1.0]:
            line.DrawLine(xmin, y, xmax/50, y)
        for i, ss in enumerate(superstrips):
            gr = graphs["gr_%s" % ss[0]]
            gr.Draw("C")
        legend.Draw()
        CMS_label()
        save(imgdir, "projective_sorted_zoom_%s" % tower)

        # Print out
        for i, ss in enumerate(superstrips):
            print '{0:7} {1:10d}  {2:5.4f}'.format(ss[0], ss[1], ss[2])
        print
        for i, ss in enumerate(superstrips):
            print '{0:7} {1:10d}  {2:5.4f}'.format(ss[0], ss[3], 0.9)

        donotdelete = [hframe]
        return donotdelete

    # Barrel 2 GeV
    bank = "patternBank_lu%s_tt27_400M.root"
    (superstrips, graphs) = bookCoverage(bank)
    d = drawCoverage(superstrips, graphs, xmax=5e7, tower="tt27")

    # Barrel 3 GeV
    #bank = "patternBank_lu%s_tt27_pt3_400M.root"
    #(superstrips, graphs) = bookCoverage(bank)
    #d = drawCoverage(superstrips, graphs, xmax=5e7, tower="tt27_pt3")