'''
import random
secert=random.randint(1,10)
temp=print("guess\n")
time=1
while  time<=3:
    guess=input()
    guess=int(guess)
    if guess==secert:
        print("right")
        break
    else:
        if guess>secert:
            print("bigger")
        else:
            print("smaller")
        time=time+1
       
print("game over")
'''
for i in range(6,-1,-1):
    print(i)
