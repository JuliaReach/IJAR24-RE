###
#   Generates the plots needed for the diagram of how consonant approximation is performed
###

using ProbabilisticReachability, ProbabilityBoundsAnalysis,
      LaTeXStrings, IntervalArithmetic
import PyPlot
using ProbabilityBoundsAnalysis.PyCall
using IntervalArithmetic: Interval
using PyPlot: figure, gca, plt

function plot_alpha_box(fe::IntervalBox, alpha::Interval, colour=[1, 1, 1], plot_alpha=0.01,
                        linewidth=2, linealpha=1)
    X = [fe.v[1].lo, fe.v[1].hi]
    Y = [fe.v[2].lo, fe.v[2].hi]
    Z = [alpha.lo, alpha.hi]

    vers = [[x, y, z] for x in X, y in Y, z in Z]

    vers = vers[:]  # 8 verts

    ## Verts of each face. 6 of them, with 4 elements each
    face_pointsX1 = [[X[1], y, z] for y in Y, z in Z][:] # 4
    face_pointsX2 = [[X[2], y, z] for y in Y, z in Z][:] # 4

    face_pointsY1 = [[x, Y[1], z] for x in X, z in Z][:] # 4
    face_pointsY2 = [[x, Y[2], z] for x in X, z in Z][:] # 4

    face_pointsZ1 = [[x, y, Z[1]] for x in X, y in Y][:] # 4
    face_pointsZ2 = [[x, y, Z[2]] for x in X, y in Y][:] # 4

    face_points = (face_pointsX1, face_pointsX2, face_pointsY1, face_pointsY2, face_pointsZ1,
                   face_pointsZ2)

    ## Edge points 12 of them, with 2 elements each

    line1 = [[X[1], X[2]], [Y[1], Y[1]], [Z[1], Z[1]]]
    line2 = [[X[1], X[2]], [Y[2], Y[2]], [Z[2], Z[2]]]
    line3 = [[X[1], X[2]], [Y[1], Y[1]], [Z[2], Z[2]]]
    line4 = [[X[1], X[2]], [Y[2], Y[2]], [Z[1], Z[1]]]
    line5 = [[X[1], X[1]], [Y[1], Y[2]], [Z[1], Z[1]]]
    line6 = [[X[2], X[2]], [Y[1], Y[2]], [Z[2], Z[2]]]
    line7 = [[X[1], X[1]], [Y[1], Y[2]], [Z[2], Z[2]]]
    line8 = [[X[2], X[2]], [Y[1], Y[2]], [Z[1], Z[1]]]
    line9 = [[X[1], X[1]], [Y[1], Y[1]], [Z[1], Z[2]]]
    line10 = [[X[2], X[2]], [Y[2], Y[2]], [Z[1], Z[2]]]
    line11 = [[X[1], X[1]], [Y[2], Y[2]], [Z[1], Z[2]]]
    line12 = [[X[2], X[2]], [Y[1], Y[1]], [Z[1], Z[2]]]

    lines = (line1, line2, line3, line4, line5, line6, line7, line8, line9, line10, line11, line12)

    ax = plt.subplot(; projection="3d")

    # plot faces
    for fp in face_points
        verts2 = ([tuple(fp[1]...); tuple(fp[3]...); tuple(fp[4]...); tuple(fp[2]...)],)

        p3c2 = PyObject(art3d.Poly3DCollection(verts2; alpha=plot_alpha))
        #p3c2 = PyObject(art3d.Poly3DCollection(verts2, alpha=0.0))
        face_color = colour
        #face_color = PyPlot.cm.RdBu(alpha.hi)
        #face_color = [alpha.hi, 0, 0]
        pycall(p3c2.set_facecolor, PyAny, face_color)
        pycall(ax.add_collection3d, PyAny, p3c2)
    end

    # plot lines

    for ln in lines
        #ax.plot3D(ln[1], ln[2], ln[3], color= PyPlot.cm.RdBu(alpha.hi), linewidth=2)
        ax.plot3D(ln[1], ln[2], ln[3]; color="black", linewidth=linewidth, alpha=linealpha)
        #ax.plot3D(ln[1], ln[2], ln[3], color= [alpha.hi, 0, 0], linewidth=2)
    end

    #ax.view_init(50,30)

    #zlim(0,1)
end

plot_dir = "$(@__DIR__)/plots"
if !isdir(plot_dir)
    mkdir(plot_dir)
end

ProbabilityBoundsAnalysis.using3D()

ft_big = 55
ft = 44
ft_ticks = 40
fig_size = (10, 10)

int_width = 4   # for line plot width
dot_width = 3   # for line plot width

Ns = 200

beta_1 = 0.2;
beta_2 = 0.6;

I_2(beta) = IntervalBox(interval(beta / 2, 1 - (beta / 2)), 2)

ProbabilityBoundsAnalysis.setSteps(Ns)
x1 = beta(interval(1, 2), 3)
x2 = gamma(interval(5, 6), 2)

### plot x1 with cuts ###

cut1 = ProbabilityBoundsAnalysis.cut(x1, beta_1 / 2 .. 1 - beta_1 / 2)
cut2 = ProbabilityBoundsAnalysis.cut(x1, beta_2 / 2 .. 1 - beta_2 / 2)

fig = figure("x1"; figsize=fig_size)
ProbabilityBoundsAnalysis.plot(x1; name="x1", fontsize=ft)
#PyPlot.xlabel(L"$X_1$", fontsize = ft)
PyPlot.xlabel("")
PyPlot.ylabel("")
#PyPlot.xticks([])
#PyPlot.yticks([])

fig = PyPlot.gcf()
#fig.set_size_inches(6.4,  4.8)

PyPlot.tight_layout()

xlimm = gca().get_xlim()
ylimm = gca().get_ylim()

xpoint = xlimm[1]
ypoint = ylimm[1]

PyPlot.plot([-0.15; cut1.hi], [1 - beta_1 / 2; 1 - beta_1 / 2], "--"; color="red",
            linewidth=dot_width)
PyPlot.plot([-0.15; cut1.lo], [beta_1 / 2; beta_1 / 2], "--"; color="red", linewidth=dot_width)
PyPlot.plot([cut1.lo; cut1.lo], [-0.15; beta_1 / 2], "--"; color="red", linewidth=dot_width)
PyPlot.plot([cut1.hi; cut1.hi], [-0.15; 1 - beta_1 / 2], "--"; color="red", linewidth=dot_width)

PyPlot.plot([-0.2; cut2.hi], [1 - beta_2 / 2; 1 - beta_2 / 2], "--"; color="blue",
            linewidth=dot_width)
PyPlot.plot([-0.2; cut2.lo], [beta_2 / 2; beta_2 / 2], "--"; color="blue", linewidth=dot_width)
PyPlot.plot([cut2.lo; cut2.lo], [-0.2; beta_2 / 2], "--"; color="blue", linewidth=dot_width)
PyPlot.plot([cut2.hi; cut2.hi], [-0.2; 1 - beta_2 / 2], "--"; color="blue", linewidth=dot_width)

PyPlot.plot([cut2.lo; cut2.hi], [-0.2; -0.2]; color="blue", linewidth=int_width)
PyPlot.plot([cut1.lo; cut1.hi], [-0.15; -0.15]; color="red", linewidth=int_width)

PyPlot.plot([-0.2; -0.2], [beta_2 / 2; 1 - beta_2 / 2]; color="blue", linewidth=int_width)
PyPlot.plot([-0.15; -0.15], [beta_1 / 2; 1 - beta_1 / 2]; color="red", linewidth=int_width)

#PyPlot.axis("off")

ax = gca()

ax.spines["right"].set_color("none")
ax.spines["top"].set_color("none")
ax.xaxis.set_ticks_position("bottom")
ax.spines["bottom"].set_position(("data", 0))
ax.yaxis.set_ticks_position("left")
ax.spines["left"].set_position(("data", 0))
ax.spines["left"].set_linewidth(2)
ax.spines["bottom"].set_linewidth(2)
ax.tick_params(; width=2)

PyPlot.xticks([0.2, 0.4, 0.6, 0.8, 1]; fontsize=ft_ticks)
PyPlot.yticks([0.2, 0.4, 0.6, 0.8, 1]; fontsize=ft_ticks)

PyPlot.title(L"$F_{X_{1}}$"; fontsize=ft_big)
PyPlot.tight_layout()

PyPlot.savefig("$plot_dir/x1_cuts.pdf")
PyPlot.savefig("$plot_dir/x1_cuts.png")

### plot x2 with cuts ###

cut1 = ProbabilityBoundsAnalysis.cut(x2, beta_1 / 2 .. 1 - beta_1 / 2)
cut2 = ProbabilityBoundsAnalysis.cut(x2, beta_2 / 2 .. 1 - beta_2 / 2)

fig = figure("x2"; figsize=fig_size)
ProbabilityBoundsAnalysis.plot(x2; name="x2", fontsize=ft)
#PyPlot.xlabel(L"$X_1$", fontsize = ft)
PyPlot.xlabel("")
PyPlot.ylabel("")
#PyPlot.xticks([])
#PyPlot.yticks([])

fig = PyPlot.gcf()
#fig.set_size_inches(6.4,  4.8)

PyPlot.tight_layout()

xlimm = gca().get_xlim()
ylimm = gca().get_ylim()

xpoint = xlimm[1]
ypoint = ylimm[1]

cut_1_corr = -right(x2) * 0.15
cut_2_corr = -right(x2) * 0.2

PyPlot.plot([cut_1_corr; cut1.hi], [1 - beta_1 / 2; 1 - beta_1 / 2], "--"; color="red",
            linewidth=dot_width)
PyPlot.plot([cut_1_corr; cut1.lo], [beta_1 / 2; beta_1 / 2], "--"; color="red", linewidth=dot_width)
PyPlot.plot([cut1.lo; cut1.lo], [-0.15; beta_1 / 2], "--"; color="red", linewidth=dot_width)
PyPlot.plot([cut1.hi; cut1.hi], [-0.15; 1 - beta_1 / 2], "--"; color="red", linewidth=dot_width)

PyPlot.plot([cut_2_corr; cut2.hi], [1 - beta_2 / 2; 1 - beta_2 / 2], "--"; color="blue",
            linewidth=dot_width)
PyPlot.plot([cut_2_corr; cut2.lo], [beta_2 / 2; beta_2 / 2], "--"; color="blue",
            linewidth=dot_width)
PyPlot.plot([cut2.lo; cut2.lo], [-0.2; beta_2 / 2], "--"; color="blue", linewidth=dot_width)
PyPlot.plot([cut2.hi; cut2.hi], [-0.2; 1 - beta_2 / 2], "--"; color="blue", linewidth=dot_width)

PyPlot.plot([cut2.lo; cut2.hi], [-0.2; -0.2]; color="blue", linewidth=int_width)
PyPlot.plot([cut1.lo; cut1.hi], [-0.15; -0.15]; color="red", linewidth=int_width)

PyPlot.plot([cut_2_corr; cut_2_corr], [beta_2 / 2; 1 - beta_2 / 2]; color="blue",
            linewidth=int_width)
PyPlot.plot([cut_1_corr; cut_1_corr], [beta_1 / 2; 1 - beta_1 / 2]; color="red",
            linewidth=int_width)

#PyPlot.axis("off")

ax = gca()

ax.spines["right"].set_color("none")
ax.spines["top"].set_color("none")
ax.xaxis.set_ticks_position("bottom")
ax.spines["bottom"].set_position(("data", 0))
ax.yaxis.set_ticks_position("left")
ax.spines["left"].set_position(("data", 0))
ax.spines["left"].set_linewidth(2)
ax.spines["bottom"].set_linewidth(2)
ax.tick_params(; width=2)

PyPlot.xticks([5, 10, 15, 20, 25, 30]; fontsize=ft_ticks)
PyPlot.yticks([0.2, 0.4, 0.6, 0.8, 1]; fontsize=ft_ticks)

PyPlot.title(L"$F_{X_{2}}$"; fontsize=ft_big)
PyPlot.tight_layout()

PyPlot.savefig("$plot_dir/x2_cuts.pdf")
PyPlot.savefig("$plot_dir/x2_cuts.png")
#PyPlot.savefig("$plot_dir/x2_no_label.png")

### make multivariate pbox ###

cor12 = -0.8;
C12 = GaussianCopula([1 cor12; cor12 1])
F12 = C12([x1, x2])

### plot X1 X2 ###

FE, masses = cons_split(F12)
alphas = [0; cumsum(masses)]
alpha_ints = interval.(alphas[1:(end - 1)], alphas[2:end])

name = "Possibility_construct"
fig = figure(name; figsize=(10, 18))
art3d = PyObject(PyPlot.art3D)
plot_alpha_box.(FE, alpha_ints, "black", Ref(0.01), Ref(1), Ref(0.5))

alpha1 = 1 - Hvolume(C12, I_2(beta_1))
alpha2 = 1 - Hvolume(C12, I_2(beta_2))

cut1_index = findlast(alphas .<= alpha1)
cut2_index = findlast(alphas .<= alpha2)

cut1 = FE[cut1_index]
cut2 = FE[cut2_index]

#alpha1_int = alpha_ints[cut1_index]
#alpha2_int = alpha_ints[cut2_index]
alpha1_int = interval(alpha1)
alpha2_int = interval(alpha2)

plot_alpha_box(cut1, alpha1_int, "red", 0.7)
plot_alpha_box(cut2, alpha2_int, "deepskyblue", 0.7)

PyPlot.zlabel(L"$\alpha$"; fontsize=ft)
PyPlot.xlabel(L"$X_1$"; fontsize=ft)
PyPlot.ylabel(L"$X_2$"; fontsize=ft)

# plots annotations

X1_top_1 = 2
X2_top_1 = -12

X1_top_2 = 1.5
X2_top_2 = -10

ax = gca()

ax.plot3D([cut1[1].hi, X1_top_1], [cut1[2].lo, cut1[2].lo], [alpha1, alpha1], "--"; color="red",
          linewidth=2)
ax.plot3D([cut1[1].hi, X1_top_1], [cut1[2].hi, cut1[2].hi], [alpha1, alpha1], "--"; color="red",
          linewidth=2)

ax.plot3D([cut1[1].hi, cut1[1].hi], [cut1[2].lo, X2_top_1], [alpha1, alpha1], "--"; color="red",
          linewidth=2)
ax.plot3D([cut1[1].lo, cut1[1].lo], [cut1[2].lo, X2_top_1], [alpha1, alpha1], "--"; color="red",
          linewidth=2)

ax.plot3D([cut2[1].hi, X1_top_2], [cut2[2].lo, cut2[2].lo], [alpha2, alpha2], "--";
          color="deepskyblue", linewidth=2)
ax.plot3D([cut2[1].hi, X1_top_2], [cut2[2].hi, cut2[2].hi], [alpha2, alpha2], "--";
          color="deepskyblue", linewidth=2)

ax.plot3D([cut2[1].hi, cut2[1].hi], [cut2[2].lo, X2_top_2], [alpha2, alpha2], "--";
          color="deepskyblue", linewidth=2)
ax.plot3D([cut2[1].lo, cut2[1].lo], [cut2[2].lo, X2_top_2], [alpha2, alpha2], "--";
          color="deepskyblue", linewidth=2)

# Fix plot params

PyPlot.xticks(; fontsize=ft_ticks)
PyPlot.yticks(; fontsize=ft_ticks)
#ax.zticks(fontsize = ft_ticks)
ax.axes.xaxis.set_ticklabels([])
ax.axes.yaxis.set_ticklabels([])

for t in ax.zaxis.get_major_ticks()
    t.label.set_fontsize(ft_ticks)
    #t.label.set_fontsize("")
end
ax.axes.zaxis.set_ticklabels([])

PyPlot.tight_layout()
ax.set_box_aspect((1, 1, 2.2))

ax.view_init(20, -42)
fig.canvas.draw()

PyPlot.savefig("$plot_dir/$name.pdf")
PyPlot.savefig("$plot_dir/$name.png")
PyPlot.close("all")
