#!/bin/sh

#  Script.sh
#  
#
#  Created by Nasuf SONMEZ on 4/28/15.
#

cat nohup.out | grep "0  event after" |  awk -F":"  '{print $2}' | awk -F " " ' {print "\t&" $2  }' > col1
cat nohup.out | grep "1  event after" |  awk -F":"  '{print $2}' | awk -F " " ' {print "\t&" $2  }' > col2
cat nohup.out | grep "2  event after" |  awk -F":"  '{print $2}' | awk -F " " ' {print "\t&" $2  }' > col3
cat nohup.out | grep "3  event after" |  awk -F":"  '{print $2}' | awk -F " " ' {print "\t&" $2  }' > col4
cat nohup.out | grep "4  event after met" |  awk -F":"  '{print $2}' | awk -F " " ' {print "\t&" $2  }' > col5
cat nohup.out | grep "5  event after" |  awk -F":"  '{print $2}' | awk -F " " ' {print "\t&" $2  }' > col5a
cat nohup.out | grep "5a event after" |  awk -F":"  '{print $2}' | awk -F " " ' {print "\t&" $2  }' > col6
#cat nohup.out | grep "6  event after" |  awk -F":"  '{print $2}' | awk -F " " ' {print "\t&" $2  }' > col7
cat nohup.out | grep "7  event after" |  awk -F":"  '{print $2}' | awk -F " " ' {print "\t&" $2  }' > col8
cat nohup.out | grep "8  event after" |  awk -F":"  '{print $2}' | awk -F " " ' {print "\t&" $2  }' > col9
cat nohup.out | grep "9  event after" |  awk -F":"  '{print $2}' | awk -F " " ' {print "\t&" $2  }' > col10
cat nohup.out | grep "10 event after"|  awk -F":"  '{print $2}' | awk -F " " ' {print "\t&" $2  }' > col11

paste col0 col1 col2 col3 col4 col5a col5 col8 col9 col10 col11 colend
