import sys
import meshio
import numpy as np

if len(sys.argv) < 3:
    print(f"Usage: {argv[0]} mesh0 mesh1 [atol]")
    exit(1)

atol = 1e-10
rtol = 1e-2

mesh0 = meshio.read(sys.argv[1])
mesh1 = meshio.read(sys.argv[2])

print(f"comparing {sys.argv[1]} and {sys.argv[2]}")

if len(sys.argv) > 3:
    atol = float(sys.argv[3])

ok = True

keys0 = mesh0.point_data.keys()
keys1 = mesh1.point_data.keys()
if keys0 != keys1:
    print(f"different data point keys: {keys0} vs {keys1}")
    ok = False
    
if mesh0.points.shape != mesh1.points.shape:
    print(f"different point array shape: {mesh0.points.shape} vs {mesh1.points.shape}")
    ok = False

if ok:
    cnt = 0
    resort0 = list(i for p, i in sorted([(tuple(mesh0.points[i, :]), i) for i in range(mesh0.points.shape[0])]))
    resort1 = list(i for p, i in sorted([(tuple(mesh1.points[i, :]), i) for i in range(mesh0.points.shape[0])]))
    
    for i in range(mesh0.points.shape[0]):
        p0 = mesh0.points[resort0[i], :]
        p1 = mesh1.points[resort1[i], :]
        if np.linalg.norm(p0-p1) > atol: #or abs(np.linalg.norm(p1)/np.linalg.norm(p0) - 1.0) > rtol:
            if ok:
                print(f"different coordinates (sorted) at position {i}: {p0} vs {p1}")
                ok = False
            cnt += 1
    if cnt > 0:
        print(f"total different coordinates (sorted): {cnt} out of {mesh0.points.shape[0]}")
    else:
        print(f"no different coordinates (sorted)")

if ok:
    print(f"keys: {keys0}")
    for k in keys0:
        cnt = 0
        resort0 = list(i for p, i in sorted([(tuple(mesh0.points[i, :]), i) for i in range(mesh0.points.shape[0])]))
        resort1 = list(i for p, i in sorted([(tuple(mesh1.points[i, :]), i) for i in range(mesh0.points.shape[0])]))
        
        diff = np.empty((mesh0.point_data[k].shape[0]), dtype=np.float64)
        if len(mesh0.point_data[k].shape) > 1 and mesh0.point_data[k].shape[1] > 1:
            print(f"key {k} is vector-valued")
            v0 = mesh0.point_data[k][resort0[0]]
            #print(f"{v0} is an example of value")
            
        for i in range(mesh0.point_data[k].shape[0]):
            v0 = mesh0.point_data[k][resort0[i]]
            v1 = mesh1.point_data[k][resort1[i]]
            if np.linalg.norm(v0-v1) > atol: #or abs(np.linalg.norm(p1)/np.linalg.norm(p0) - 1.0) > rtol:
                if ok:
                    print(f"different {k}-values (sorted) at position {i}: {v0} vs {v1}")
                    ok = False
                cnt += 1
                diff[i] = np.linalg.norm(v0-v1)
            else:
                diff[i] = 0
        
        if cnt > 0:
            print(f"total different {k}-values (sorted): {cnt} out of {mesh0.point_data[k].shape[0]}")
            if abs(np.var(diff)) < atol:
                print(f"variance of difference is negligible, mean is {np.mean(diff)}")
        else:
            print(f"no different {k}-values (sorted)")
