using ProbabilisticReachability, ProbabilityBoundsAnalysis,
      LaTeXStrings, IntervalArithmetic
import PyPlot
using ProbabilityBoundsAnalysis.PyCall
using IntervalArithmetic: Interval
using PyPlot: figure, gca, plt

function islesseqhi(box::IntervalBox, p::Vector)
    return (box.v[1].hi <= p[1]) && (box.v[2].hi <= p[2])
end

function islesseqlo(box::IntervalBox, p::Vector)
    return (box.v[1].lo <= p[1]) && (box.v[2].lo <= p[2])
end

function plot_splits_2D(F)
    fe, masses = split(F)
    masses = masses[:]

    cdf_bounds = Vector{Interval{Float64}}(undef, length(fe))

    for (i, f_e) in enumerate(fe)
        p1 = [f_e.v[1].lo, f_e.v[2].lo]
        p2 = [f_e.v[1].hi, f_e.v[2].hi]

        these1 = islesseqhi.(fe, Ref(p1))
        these2 = islesseqhi.(fe, Ref(p2))

        pls = sum(masses[these2])
        bel = sum(masses[these1])

        cdf_bounds[i] = interval(bel, pls)
    end

    return plot_cdf_box.(fe, cdf_bounds)
end

function plot_cdf_box(fe::IntervalBox, cdf_bound::Interval)
    X = [fe.v[1].lo, fe.v[1].hi]
    Y = [fe.v[2].lo, fe.v[2].hi]
    Z = [cdf_bound.lo, cdf_bound.hi]

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

        p3c2 = PyObject(art3d.Poly3DCollection(verts2; alpha=0.05))
        face_color = [1, 0, 0]
        pycall(p3c2.set_facecolor, PyAny, face_color)
        pycall(ax.add_collection3d, PyAny, p3c2)
    end

    # plot lines

    for ln in lines
        ax.plot3D(ln[1], ln[2], ln[3]; color="b", linewidth=2)
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

ProbabilityBoundsAnalysis.setSteps(5)
x1 = beta(interval(1, 2), 3)
ProbabilityBoundsAnalysis.setSteps(3)
x2 = gamma(interval(5, 6), 2)
ProbabilityBoundsAnalysis.setSteps(10)
x3 = normal(interval(0, 0.5), interval(2, 3))

cor12 = -0.8;
cor13 = 0.8;
cor23 = -0.8;

C12 = GaussianCopula([1 cor12; cor12 1])
C13 = GaussianCopula([1 cor13; cor13 1])
C23 = GaussianCopula([1 cor23; cor23 1])

F12 = C12([x1, x2])
F13 = C12([x1, x3])
F23 = C12([x2, x3])

name = "biv_x1x2_split"
fig = figure(name; figsize=fig_size)
art3d = PyObject(PyPlot.art3D)
plot_splits_2D(F12)

PyPlot.zlabel(L"$F$"; fontsize=ft)
PyPlot.xlabel(L"$X_1$"; fontsize=ft)
PyPlot.ylabel(L"$X_2$"; fontsize=ft)

ax = gca()
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

ax.view_init(30, -60)
fig.canvas.draw()
PyPlot.title(L"$F_{X_1,X_2}$"; fontsize=ft_big)
PyPlot.tight_layout()
PyPlot.savefig("$plot_dir/$name.pdf")
PyPlot.savefig("$plot_dir/$name.png")
PyPlot.close("all")

name = "biv_x1x3_split"
fig = figure(name; figsize=(10, 10))
art3d = PyObject(PyPlot.art3D)
plot_splits_2D(F13)

PyPlot.zlabel(L"$F$"; fontsize=ft)
PyPlot.xlabel(L"$X_1$"; fontsize=ft)
PyPlot.ylabel(L"$X_3$"; fontsize=ft)

ax = gca()
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

ax.view_init(25, -60)
fig.canvas.draw()
PyPlot.title(L"$F_{X_1, X_3}$"; fontsize=ft_big)
PyPlot.tight_layout()

PyPlot.savefig("$plot_dir/$name.pdf")
PyPlot.savefig("$plot_dir/$name.png")
PyPlot.close("all")

name = "biv_x2x3_split"
fig = figure(name; figsize=fig_size)
art3d = PyObject(PyPlot.art3D)
plot_splits_2D(F23)

PyPlot.zlabel(L"$F$"; fontsize=ft)
PyPlot.xlabel(L"$X_2$"; fontsize=ft)
PyPlot.ylabel(L"$X_3$"; fontsize=ft)

ax = gca()
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

ax.view_init(30, -60)
fig.canvas.draw()
PyPlot.title(L"$F_{X_2,X_3}$"; fontsize=ft_big)
PyPlot.tight_layout()
PyPlot.savefig("$plot_dir/$name.pdf")
PyPlot.savefig("$plot_dir/$name.png")
PyPlot.close("all")

#######################

name = "biv_x1x2_2_split"
fig = figure(name; figsize=fig_size)
art3d = PyObject(PyPlot.art3D)
plot_splits_2D(F12)

PyPlot.zlabel(L"$F$"; fontsize=ft)
PyPlot.xlabel(L"$X_1$"; fontsize=ft)
PyPlot.ylabel(L"$X_2$"; fontsize=ft)

ax = gca()
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

ax.view_init(30, 180 - 30)
fig.canvas.draw()
PyPlot.title(L"$F_{X_1, X_2}$"; fontsize=ft_big)
PyPlot.tight_layout()
PyPlot.savefig("$plot_dir/$name.pdf")
PyPlot.savefig("$plot_dir/$name.png")
PyPlot.close("all")

name = "biv_x1x3_2_split"
fig = figure(name; figsize=fig_size)

art3d = PyObject(PyPlot.art3D)
plot_splits_2D(F13)

PyPlot.zlabel(L"$F$"; fontsize=ft)
PyPlot.xlabel(L"$X_1$"; fontsize=ft)
PyPlot.ylabel(L"$X_3$"; fontsize=ft)

ax = gca()
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

ax.view_init(30, 180 - 30)
fig.canvas.draw()
PyPlot.title(L"$F_{X_1,X_3}$"; fontsize=ft_big)
PyPlot.tight_layout()
PyPlot.savefig("$plot_dir/$name.pdf")
PyPlot.savefig("$plot_dir/$name.png")
PyPlot.close("all")

name = "biv_x2x3_2_split"
fig = figure(name; figsize=fig_size)
art3d = PyObject(PyPlot.art3D)
plot_splits_2D(F23)

PyPlot.zlabel(L"$F$"; fontsize=ft)

PyPlot.xlabel(L"$X_2$"; fontsize=ft)
PyPlot.ylabel(L"$X_3$"; fontsize=ft)

ax = gca()
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

ax.view_init(30, 180 - 30)
fig.canvas.draw()
PyPlot.title(L"$F_{X_2,X_3}$"; fontsize=ft_big)
PyPlot.tight_layout()
PyPlot.savefig("$plot_dir/$name.pdf")
PyPlot.savefig("$plot_dir/$name.png")
PyPlot.close("all")

#=
    feLazy = convert.(Hyperrectangle, fe)
    for (i, f_e) in enumerate(feLazy)

        indbel = [issubset(feL, f_e) for feL in feLazy];
        indpls = [!isdisjoint(feL, f_e) for feL in feLazy];

        bel = sum(masses[indbel]);
        pls = sum(masses[indpls]);

        cdf_bounds[i] = interval(bel, pls)

    end
    =#
#=  Looked good, but doesn't work

Cop = F.copula
n1 = F.marginals[1].n
n2 = F.marginals[2].n

i1 = range(0, 1, length = n1+1)
i2 = range(0, 1, length = n2+1)

islo = [[ii1, ii2] for ii1 in i1[1:end-1], ii2 in i2[1:end-1]]
ishi = [[ii1, ii2] for ii1 in i1[2:end], ii2 in i2[2:end]]

cdfslo = [cdf(Cop,i) for i in islo]
cdfshi = [cdf(Cop,i) for i in ishi]

cdf_bounds = interval.(cdfslo[:], cdfshi[:])
=#
