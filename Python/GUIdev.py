#%%
#---------Imports
from numpy import arange, sin, pi
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import tkinter as Tk
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
#---------End of imports

fig = plt.Figure()

x = np.arange(0, 2*np.pi, 0.01)        # x-array

def animate(i):
    line.set_ydata(np.sin(x+i/10.0))  # update the data
    return line,

root = Tk.Tk()

label = Tk.Label(root,text="SHM Simulation").grid(column=0, row=0)

canvas = FigureCanvasTkAgg(fig, master=root)
canvas.get_tk_widget().grid(column=0,row=1)

ax = fig.add_subplot(111)
line, = ax.plot(x, np.sin(x))
ani = animation.FuncAnimation(fig, animate, np.arange(1, 200), interval=25, blit=False)

Tk.mainloop()
#%%
#---------Imports
from numpy import arange, sin, pi
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import tkinter as Tk
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
#---------End of imports

fig = plt.Figure()

x = np.arange(0, 2*np.pi, 0.01)        # x-array

def animate(i):
    line.set_ydata(np.sin(x+i/10.0))  # update the data
    return line,

root = Tk.Tk()

label = Tk.Label(root,text="SHM Simulation").grid(column=0, row=0)

canvas = FigureCanvasTkAgg(fig, master=root)
canvas.get_tk_widget().grid(column=0,row=1)
#canvas = FigureCanvasTkAgg(fig, root)
canvas.draw()
canvas.get_tk_widget().pack(fill=Tk.BOTH, expand=True)

ax = fig.add_subplot(111)
line, = ax.plot(x, np.sin(x))
ani = animation.FuncAnimation(fig, animate, np.arange(1, 200), interval=25, blit=False)

Tk.mainloop()

#%%
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg 
import tkinter as tk
import matplotlib.pyplot as plt
 
root = tk.Tk()
fig, ax = plt.subplots()

t = np.arange(0, 3, .01)
line, = ax.plot(t, 2 * np.sin(2 * np.pi * t))


canvas = FigureCanvasTkAgg(fig, root)
canvas.draw()
#canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
canvas.get_tk_widget().pack(side=tk.LEFT)

canvas.get_tk_widget().pack(side=tk.LEFT)




 
root.mainloop()

#%%

#---------Imports
from numpy import arange, sin, pi
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import tkinter as tk
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
#---------End of imports



fig = plt.Figure()

x = np.arange(0, 2*np.pi, 0.01)        # x-array

def animate(i):
    line.set_ydata(np.sin(x+i/10.0))  # update the data
    return line,


root = tk.Tk()
canvas = tk.Canvas(root,width=600,height=300)
canvas.grid(columnspan=3)


ax = fig.add_subplot(111)
line, = ax.plot(x, np.sin(x))
fig1 = animation.FuncAnimation(fig, animate, np.arange(1, 200), interval=25, blit=False)


fig_label = tk.Label(fig=fig)
fig_label.image=fig1
fig_label.grid(column=1,row=0)
root.mainloop()

#%%
import tkinter as tk
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


#data input
data1 = {'country': ['A', 'B', 'C', 'D', 'E'],
         'gdp_per_capita': [45000, 42000, 52000, 49000, 47000]
         }
df1 = pd.DataFrame(data1)

data2 = {'year': [1920, 1930, 1940, 1950, 1960, 1970, 1980, 1990, 2000, 2010],
         'unemployment_rate': [9.8, 12, 8, 7.2, 6.9, 7, 6.5, 6.2, 5.5, 6.3]
         }  
df2 = pd.DataFrame(data2)

data3 = {'interest_rate': [5, 5.5, 6, 5.5, 5.25, 6.5, 7, 8, 7.5, 8.5],
         'index_price': [1500, 1520, 1525, 1523, 1515, 1540, 1545, 1560, 1555, 1565]
         }
df3 = pd.DataFrame(data3)


#tkinter startloop
root = tk.Tk()




#Figure 1
figure1 = plt.Figure(figsize=(6, 6), dpi=100)
ax1 = figure1.add_subplot(111)

bar1 = FigureCanvasTkAgg(figure1, root)
bar1.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH)

df1 = df1[['country', 'gdp_per_capita']].groupby('country').sum()
df1.plot(kind='bar', legend=True, ax=ax1)
ax1.set_title('Country Vs. GDP Per Capita')

#Figure 2
figure2 = plt.Figure(figsize=(5, 2), dpi=100)
ax2 = figure2.add_subplot(111)

line2 = FigureCanvasTkAgg(figure2, root)
line2.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH)

df2 = df2[['year', 'unemployment_rate']].groupby('year').sum()
df2.plot(kind='line', legend=True, ax=ax2, color='r', marker='o', fontsize=10)
ax2.set_title('Year Vs. Unemployment Rate')

#Figure 3
figure3 = plt.Figure(figsize=(5, 2), dpi=100)
ax3 = figure3.add_subplot(111)
ax3.scatter(df3['interest_rate'], df3['index_price'], color='g')

scatter3 = FigureCanvasTkAgg(figure3, root)
scatter3.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH)

ax3.legend(['index_price'])
ax3.set_xlabel('Interest Rate')
ax3.set_title('Interest Rate Vs. Index Price')

#Figure 4
figure4 = plt.Figure(figsize=(5, 3), dpi=100)
ax4 = figure4.add_subplot(111)
ax4.scatter(df3['interest_rate'], df3['index_price'], color='g')

scatter4 = FigureCanvasTkAgg(figure4, root)
scatter4.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH)

ax4.legend(['index_price'])
ax4.set_xlabel('Interest Rate')
ax4.set_title('Interest Rate Vs. Index Price')


root.mainloop()

#%%

import tkinter as tk
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


data3 = {'interest_rate': [5, 5.5, 6, 5.5, 5.25, 6.5, 7, 8, 7.5, 8.5],
         'index_price': [1500, 1520, 1525, 1523, 1515, 1540, 1545, 1560, 1555, 1565]
         }
df3 = pd.DataFrame(data3)


root = tk.Tk()
root.title('ROXS GUI')

frame=tk.Frame(root)
frame.grid(column=0, row=0)

#current graph
l0 = tk.LabelFrame(frame,text="Current Plot")
l0.grid(row=0,column=0)

l0label=tk.Label(l0,text='Current Plot1')
l0label.grid(row=0,column=0)

#stored graphs
l1 = tk.Frame(frame)
l1.grid(row=0,column=1)


l1label=tk.Label(l1,text='Stored Plots1')
l1label.grid(row=1,column=0)

l10 = tk.LabelFrame(l1,text='Spiral Graphs').grid(column=0, row=0)

l10label=tk.Label(l10,text='plotplot')
#l10label.grid(row=0,column=0)

#l11 = tk.Frame(l1).grid(column=0, row=1)

#l12 = tk.Frame(l1).grid(column=0, row=2)


#fig=tk.Entry(l10).grid(row=1, column=0)

# #Spiral Figures
# figure3 = plt.Figure(figsize=(5, 2), dpi=100)
# ax3 = figure3.add_subplot(111)
# ax3.scatter(df3['interest_rate'], df3['index_price'], color='g')

# scatter3 = FigureCanvasTkAgg(figure3, l10)
# scatter3.get_tk_widget().grid(column=0, row=0)

# scatter4 = FigureCanvasTkAgg(figure3, l1)
# scatter4.get_tk_widget().grid(column=1, row=0)

# scatter5 = FigureCanvasTkAgg(figure3, l1)
# scatter5.get_tk_widget().grid(column=2, row=0)

# #Linear Figures
# figure3 = plt.Figure(figsize=(5, 2), dpi=100)
# ax3 = figure3.add_subplot(111)
# ax3.scatter(df3['interest_rate'], df3['index_price'], color='g')

# scatter3 = FigureCanvasTkAgg(figure3, l1)
# scatter3.get_tk_widget().grid(column=0, row=0)

# scatter4 = FigureCanvasTkAgg(figure3, l1)
# scatter4.get_tk_widget().grid(column=1, row=0)

# scatter5 = FigureCanvasTkAgg(figure3, l1)
# scatter5.get_tk_widget().grid(column=2, row=0)

# scatter3 = FigureCanvasTkAgg(figure3, l11)
# scatter3.get_tk_widget().grid(row=1,column=0)

# scatter4 = FigureCanvasTkAgg(figure3, l11)
# scatter4.get_tk_widget().grid(row=1,column=1)

# scatter5 = FigureCanvasTkAgg(figure3, l11)
# scatter5.get_tk_widget().grid(row=1,column=2)

# l11 = tk.Label(l1,text="Linear Graphs").grid(column=0, row=1)
# l12 = tk.Label(l1,text="Push Pull Graphs").grid(column=0, row=2)



root.mainloop()

#%%
import tkinter
from tkinter import ttk
from tkinter import messagebox

def enter_data():
    accepted = accept_var.get()
    
    if accepted=="Accepted":
        # User info
        firstname = first_name_entry.get()
        lastname = last_name_entry.get()
        
        if firstname and lastname:
            title = title_combobox.get()
            age = age_spinbox.get()
            nationality = nationality_combobox.get()
            
            # Course info
            registration_status = reg_status_var.get()
            numcourses = numcourses_spinbox.get()
            numsemesters = numsemesters_spinbox.get()
            
            print("First name: ", firstname, "Last name: ", lastname)
            print("Title: ", title, "Age: ", age, "Nationality: ", nationality)
            print("# Courses: ", numcourses, "# Semesters: ", numsemesters)
            print("Registration status", registration_status)
            print("------------------------------------------")
        else:
            tkinter.messagebox.showwarning(title="Error", message="First name and last name are required.")
    else:
        tkinter.messagebox.showwarning(title= "Error", message="You have not accepted the terms")

window = tkinter.Tk()
window.title("Data Entry Form")

ww=tkinter.Frame(window)
ww.pack()



frame = tkinter.Frame(ww)
frame.pack()

frame2 = tkinter.Frame(ww)
frame2.pack()
first = tkinter.Label(user_info_frame, text="First Name")
first.grid(row=0, column=0)





# Saving User Info
user_info_frame =tkinter.LabelFrame(frame, text="User Information")
user_info_frame.grid(row= 0, column=0, padx=20, pady=10)

first_name_label = tkinter.Label(user_info_frame, text="First Name")
first_name_label.grid(row=0, column=0)
last_name_label = tkinter.Label(user_info_frame, text="Last Name")
last_name_label.grid(row=0, column=1)

first_name_entry = tkinter.Entry(user_info_frame).grid(row=1, column=0)
last_name_entry = tkinter.Entry(user_info_frame)
#first_name_entry.grid(row=1, column=0)
last_name_entry.grid(row=1, column=1)

title_label = tkinter.Label(user_info_frame, text="Title")
title_combobox = ttk.Combobox(user_info_frame, values=["", "Mr.", "Ms.", "Dr."])
title_label.grid(row=0, column=2)
title_combobox.grid(row=1, column=2)

age_label = tkinter.Label(user_info_frame, text="Age")
age_spinbox = tkinter.Spinbox(user_info_frame, from_=18, to=110)
age_label.grid(row=2, column=0)
age_spinbox.grid(row=3, column=0)

nationality_label = tkinter.Label(user_info_frame, text="Nationality")
nationality_combobox = ttk.Combobox(user_info_frame, values=["Africa", "Antarctica", "Asia", "Europe", "North America", "Oceania", "South America"])
nationality_label.grid(row=2, column=1)
nationality_combobox.grid(row=3, column=1)

for widget in user_info_frame.winfo_children():
    widget.grid_configure(padx=10, pady=5)

# Saving Course Info
courses_frame = tkinter.LabelFrame(frame)
courses_frame.grid(row=1, column=0, sticky="news", padx=20, pady=10)

registered_label = tkinter.Label(courses_frame, text="Registration Status")

reg_status_var = tkinter.StringVar(value="Not Registered")
registered_check = tkinter.Checkbutton(courses_frame, text="Currently Registered",
                                       variable=reg_status_var, onvalue="Registered", offvalue="Not registered")

registered_label.grid(row=0, column=0)
registered_check.grid(row=1, column=0)

numcourses_label = tkinter.Label(courses_frame, text= "# Completed Courses")
numcourses_spinbox = tkinter.Spinbox(courses_frame, from_=0, to='infinity')
numcourses_label.grid(row=0, column=1)
numcourses_spinbox.grid(row=1, column=1)

numsemesters_label = tkinter.Label(courses_frame, text="# Semesters")
numsemesters_spinbox = tkinter.Spinbox(courses_frame, from_=0, to="infinity")
numsemesters_label.grid(row=0, column=2)
numsemesters_spinbox.grid(row=1, column=2)

for widget in courses_frame.winfo_children():
    widget.grid_configure(padx=10, pady=5)

# Accept terms
terms_frame = tkinter.LabelFrame(frame, text="Terms & Conditions")
terms_frame.grid(row=2, column=0, sticky="news", padx=20, pady=10)

accept_var = tkinter.StringVar(value="Not Accepted")
terms_check = tkinter.Checkbutton(terms_frame, text= "I accept the terms and conditions.",
                                  variable=accept_var, onvalue="Accepted", offvalue="Not Accepted")
terms_check.grid(row=0, column=0)

# Button
button = tkinter.Button(frame, text="Enter data", command= enter_data)
button.grid(row=3, column=0, sticky="news", padx=20, pady=10)
 
window.mainloop()