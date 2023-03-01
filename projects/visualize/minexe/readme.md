# 💻🚔

The script ```scripter.m``` is compiled with,

in 🚀 & 💩
```matlab
>> mcc -d exe -m scripter.m
```

which saves the executable (along with other 🧻) in ```exe/```.

## 🏃

You are in ```exe/``` now.

in 💩
* by clicking ```scripter.exe```, the output is ```../pics/somepic.png```.
* you can also run it in *Power💩* with ```.\scripter.exe```.

in 🚀
```
$ ./run_scripter.sh ../../../../../matlabruntime/theruntimeishere/R2022b/
```
## 🧠

if you want the ```.exe``` to read ```.txt``` files,

```matlab
pwd_=pwd;
fid = fopen(fullfile(pwd_,'file.txt'));
... stuff ...
fclose(fid);

filemat = load(fullfile(pwd_,'file.txt'));
```

## 🎨

[![](../pics/somepic.png)](./)
