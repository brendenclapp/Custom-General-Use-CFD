import tkinter as tk
import numpy as np


def GUI(geo, Faces):

    Nx, Ny = geo.Nx, geo.Ny
    cell = 50
    pad = 20

    # Cell codes:
    # True  = fluid
    # False = wall
    cells =  Faces.cells

    # Face codes:
    # 0 = fluid
    # 1 = wall
    # 2 = inlet
    # 3 = outlet
    NS_faces = Faces.NS_faces
    WE_faces = Faces.WE_faces

    root = tk.Tk()
    root.title("CFD Geometry Editor")

    canvas = tk.Canvas(
        root,
        width=Nx * cell + 2 * pad,
        height=Ny * cell + 2 * pad,
        bg="white"
    )
    canvas.pack()

    NS_draw = {}
    WE_draw = {}

    def face_color(state):
        return {
            0: "#B0B0B0",  # fluid
            1: "black",   # wall
            2: "blue",    # inlet
            3: "green"    # outlet
        }[state]

    def style_face(item, state, hover=False):
        if hover:
            canvas.itemconfig(item, fill="orange", width=5)
        else:
            canvas.itemconfig(item, fill=face_color(state), width=4)

    def update_NS(i, j):
        style_face(NS_draw[(i, j)], NS_faces[i, j])

    def update_WE(i, j):
        style_face(WE_draw[(i, j)], WE_faces[i, j])

    def set_cell_faces(i, j, state):
        # bottom face
        NS_faces[i, j] = state
        update_NS(i, j)

        # top face
        NS_faces[i, j + 1] = state
        update_NS(i, j + 1)

        # left face
        WE_faces[i, j] = state
        update_WE(i, j)

        # right face
        WE_faces[i + 1, j] = state
        update_WE(i + 1, j)

    def toggle_cell(i, j, rect):
        cells[i, j] = not cells[i, j]

        if cells[i, j]:
            canvas.itemconfig(rect, fill="white")
            set_cell_faces(i, j, 0)   # fluid, overwrites inlet/outlet
        else:
            canvas.itemconfig(rect, fill="#303030")
            set_cell_faces(i, j, 1)   # wall, overwrites inlet/outlet

    def cycle_face(arr, draw_dict, i, j):
        # Face clicks only cycle inlet/outlet.
        # Any current state becomes inlet first.
        if arr[i, j] == 2:
            arr[i, j] = 3
        else:
            arr[i, j] = 2

        style_face(draw_dict[(i, j)], arr[i, j])

    def make_NS_face(x1, y1, x2, y2, i, j):
        draw = canvas.create_line(
            x1, y1, x2, y2,
            fill=face_color(NS_faces[i, j]),
            width=4
        )
        NS_draw[(i, j)] = draw

        hitbox = canvas.create_line(
            x1, y1, x2, y2,
            fill="",
            width=20
        )

        canvas.tag_bind(hitbox, "<Button-1>", lambda e, i=i, j=j: cycle_face(NS_faces, NS_draw, i, j))
        canvas.tag_bind(hitbox, "<Enter>", lambda e, i=i, j=j: style_face(NS_draw[(i, j)], NS_faces[i, j], True))
        canvas.tag_bind(hitbox, "<Leave>", lambda e, i=i, j=j: style_face(NS_draw[(i, j)], NS_faces[i, j]))

    def make_WE_face(x1, y1, x2, y2, i, j):
        draw = canvas.create_line(
            x1, y1, x2, y2,
            fill=face_color(WE_faces[i, j]),
            width=4
        )
        WE_draw[(i, j)] = draw

        hitbox = canvas.create_line(
            x1, y1, x2, y2,
            fill="",
            width=20
        )

        canvas.tag_bind(hitbox, "<Button-1>", lambda e, i=i, j=j: cycle_face(WE_faces, WE_draw, i, j))
        canvas.tag_bind(hitbox, "<Enter>", lambda e, i=i, j=j: style_face(WE_draw[(i, j)], WE_faces[i, j], True))
        canvas.tag_bind(hitbox, "<Leave>", lambda e, i=i, j=j: style_face(WE_draw[(i, j)], WE_faces[i, j]))

    # cells
    # Now cells[0,0] appears at the bottom-left of the GUI
    for i in range(Nx):
        for j in range(Ny):
            x1 = pad + i * cell
            x2 = pad + (i + 1) * cell

            y1 = pad + (Ny - 1 - j) * cell
            y2 = pad + (Ny - j) * cell

            rect = canvas.create_rectangle(
                x1, y1, x2, y2,
                fill="white",
                outline="#A0A0A0",
                width=1
            )

            canvas.tag_bind(
                rect,
                "<Button-1>",
                lambda e, i=i, j=j, r=rect: toggle_cell(i, j, r)
            )

    # NS faces: horizontal faces
    # NS_faces[i,0] is the bottom face row
    # NS_faces[i,Ny] is the top face row
    for i in range(Nx):
        for j in range(Ny + 1):
            x1 = pad + i * cell
            x2 = pad + (i + 1) * cell

            y = pad + (Ny - j) * cell

            make_NS_face(x1, y, x2, y, i, j)

    # WE faces: vertical faces
    # WE_faces[0,j] is the left face column
    # WE_faces[Nx,j] is the right face column
    for i in range(Nx + 1):
        for j in range(Ny):
            x = pad + i * cell

            y1 = pad + (Ny - 1 - j) * cell
            y2 = pad + (Ny - j) * cell

            make_WE_face(x, y1, x, y2, i, j)

    def print_all():
        print("Cells:")
        print("True = fluid, False = wall")
        print(np.flipud(cells.T))

        print("\nNS_faces / horizontal faces:")
        print("0 = fluid, 1 = wall, 2 = inlet, 3 = outlet")
        print(np.flipud(NS_faces.T))

        print("\nWE_faces / vertical faces:")
        print("0 = fluid, 1 = wall, 2 = inlet, 3 = outlet")
        print(np.flipud(WE_faces.T))

    tk.Button(root, text="Print", command=print_all).pack()

    root.mainloop()

    return
