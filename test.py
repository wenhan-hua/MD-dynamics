v=trajectory_v[total_steps]
vx=[]
for t in range(len(v)):
    vx.append(v[t][0])
plt.figure()
plt.title('Distribution of the velocity')
plt.hist(vx)


plt.show()

C=c(U,T)
print(C)
