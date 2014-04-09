import sys, os

total_line_count = 0

counts = dict()

def scan_dir(dir):
  global total_line_count
  global counts
  for root, dirs, files in os.walk(dir):
      if root.__contains__('packages'): continue
      for f in files:
          # this is just raw data, skip it
          tk = os.path.splitext(f)
          ext = tk[len(tk)-1]
          if ((ext == '.dart') or
              (ext == '.html') or
              (ext == '.css') or
              (ext == '.php') or
              (ext == '.js')):
              fname = root + '/' + f
              num_lines = sum(1 for line in open(fname))
              total_line_count += num_lines
              if not counts.has_key(num_lines):
                  counts[num_lines] = [fname]
              else:
                  counts[num_lines] += [fname]

scan_dir('lib')

print('----------------------')
lc = counts.keys()
lc.sort()
for k in lc:
    for f in counts[k]:
        print k, ':', f

print '------------------'
print 'Total Lines:', total_line_count

#num_lines = sum(1 for line in open('Milo.py'))
#print num_lines
