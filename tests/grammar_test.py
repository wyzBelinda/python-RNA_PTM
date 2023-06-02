import re

a='abcdefg'
print(a[1:3])

my_dict = {'apple': 1, 'banana': 2, 'orange': 3}

if 'apple' in my_dict:
    print('The key "apple" exists in the dictionary')
else:
    print('The key "apple" does not exist in the dictionary')

sq="ADSEFRTG"
print(sq[:-1][::-1]+sq[-1])

print(1/2)
print(6/3/2)
print(6-2+2*4)

print("asd" in "asdfg")

print("a"*10)


a="02_OH-AGUC-OH.381.-1.CID_20"
b="GA"
print(b[::-1])
print(b[::-1] in a)


s = "a_11^21-C^11"
result = re.findall(r'(\w)', s)[0]
result1 = re.findall(r'_(\d+)', s)[0]
result2 = re.findall(r'\^(\d+)', s)[0]

print(result)
print(result1)
print(result2)
