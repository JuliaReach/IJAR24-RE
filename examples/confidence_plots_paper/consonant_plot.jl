###
#   Generates the 3 two-dimensional confidence distributions from the multivariate p-box
###

using ProbabilisticReachability, ProbabilityBoundsAnalysis,
      LaTeXStrings, IntervalArithmetic
import PyPlot
using ProbabilityBoundsAnalysis.PyCall
using IntervalArithmetic: Interval
using PyPlot: figure, gca, plt

function plot_alpha_box(fe::IntervalBox, alpha::Interval)
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

        p3c2 = PyObject(art3d.Poly3DCollection(verts2; alpha=0.1))
        #p3c2 = PyObject(art3d.Poly3DCollection(verts2, alpha=0.0))
        #face_color = [1, 0, 0]
        face_color = PyPlot.cm.RdBu(alpha.hi)
        #face_color = [alpha.hi, 0, 0]
        pycall(p3c2.set_facecolor, PyAny, face_color)
        pycall(ax.add_collection3d, PyAny, p3c2)
    end

    # plot lines

    for ln in lines
        #ax.plot3D(ln[1], ln[2], ln[3], color= PyPlot.cm.RdBu(alpha.hi), linewidth=2)
        ax.plot3D(ln[1], ln[2], ln[3]; color="black", linewidth=2)
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

Ns = 100

ProbabilityBoundsAnalysis.setSteps(Ns)
x1 = beta(interval(1, 2), 3)
x2 = gamma(interval(5, 6), 2)
x3 = normal(interval(0, 0.5), interval(2, 3))

cor12 = -0.8;
cor13 = 0.8;
cor23 = -0.8;

C12 = GaussianCopula([1 cor12; cor12 1])
C13 = GaussianCopula([1 cor13; cor13 1])
C23 = GaussianCopula([1 cor23; cor23 1])

F12 = C12([x1, x2])
F23 = C12([x2, x3])
F31 = C12([x3, x1])

### plot X1 X2 ###

FE, masses = cons_split(F12)
alphas = [0; cumsum(masses)]
alpha_ints = interval.(alphas[1:(end - 1)], alphas[2:end])

name = "Possibility_1"
fig = figure(name; figsize=fig_size)
art3d = PyObject(PyPlot.art3D)
plot_alpha_box.(FE, alpha_ints)

PyPlot.zlabel(L"$\alpha$"; fontsize=ft)
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

PyPlot.tight_layout()
PyPlot.savefig("$plot_dir/$name.pdf")
PyPlot.savefig("$plot_dir/$name.png")
PyPlot.close("all")

### plot X2 X3 ###

FE, masses = cons_split(F23)
alphas = [0; cumsum(masses)]
alpha_ints = interval.(alphas[1:(end - 1)], alphas[2:end])

name = "Possibility_2"
fig = figure(name; figsize=fig_size)
art3d = PyObject(PyPlot.art3D)
plot_alpha_box.(FE, alpha_ints)

PyPlot.zlabel(L"$\alpha$"; fontsize=ft)
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

PyPlot.tight_layout()
PyPlot.savefig("$plot_dir/$name.pdf")
PyPlot.savefig("$plot_dir/$name.png")
PyPlot.close("all")

### plot X1 X3 ###

FE, masses = cons_split(F31)
alphas = [0; cumsum(masses)]
alpha_ints = interval.(alphas[1:(end - 1)], alphas[2:end])

name = "Possibility_3"
fig = figure(name; figsize=fig_size)
art3d = PyObject(PyPlot.art3D)
plot_alpha_box.(FE, alpha_ints)

PyPlot.zlabel(L"$\alpha$"; fontsize=ft)
PyPlot.xlabel(L"$X_3$"; fontsize=ft)
PyPlot.ylabel(L"$X_1$"; fontsize=ft)

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

PyPlot.tight_layout()
PyPlot.savefig("$plot_dir/$name.pdf")
PyPlot.savefig("$plot_dir/$name.png")
PyPlot.close("all")
