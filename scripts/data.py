import numpy as np
import pandas as pd
import re
import datetime as dtm

p = [
    "上海2022年5月15日，新增本土新冠肺炎确诊病例69例 新增本土无症状感染者869例 新增境外输入性新冠肺炎确诊病例1例 无新增境外输入性无症状感染者",
    "【最新】5月14日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年5月14日，新增本土新冠肺炎确诊病例166例 新增本土无症状感染者1203例 新增境外输入性新冠肺炎确诊病例2例 新增境外输入性无症状感染者2例",
    "5月13日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年5月13日，新增本土新冠肺炎确诊病例194例 新增本土无症状感染者1487例 新增境外输入性新冠肺炎确诊病例1例 无新增境外输入性无症状感染者",
    "【最新】5月12日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年5月12日，新增本土新冠肺炎确诊病例227例 新增本土无症状感染者1869例 无新增境外输入性新冠肺炎确诊病例 新增境外输入性无症状感染者2例",
    "5月11日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年5月11日，新增本土新冠肺炎确诊病例144例 新增本土无症状感染者1305例 无新增境外输入性新冠肺炎确诊病例 新增境外输入性无症状感染者1例",
    "5月10日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年5月10日，新增本土新冠肺炎确诊病例228例 新增本土无症状感染者1259例 新增境外输入性新冠肺炎确诊病例1例 新增境外输入性无症状感染者1例",
    "5月9日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年5月9日，新增本土新冠肺炎确诊病例234例 新增本土无症状感染者2780例 无新增境外输入性新冠肺炎确诊病例 新增境外输入性无症状感染者3例",
    "5月8日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年5月8日，新增本土新冠肺炎确诊病例322例 新增本土无症状感染者3625例 新增境外输入性新冠肺炎确诊病例2例 新增境外输入性无症状感染者1例",
    "5月7日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年5月7日，新增本土新冠肺炎确诊病例215例 新增本土无症状感染者3760例 新增境外输入性新冠肺炎确诊病例2例 无新增境外输入性无症状感染者",
    "【最新】5月6日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年5月6日，新增本土新冠肺炎确诊病例253例 新增本土无症状感染者3961例 无新增境外输入性新冠肺炎确诊病例 无新增境外输入性无症状感染者",
    "5月5日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年5月5日，新增本土新冠肺炎确诊病例245例 新增本土无症状感染者4024例 无新增境外输入性新冠肺炎确诊病例 新增境外输入性无症状感染者2例",
    "5月4日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年5月4日，新增本土新冠肺炎确诊病例261例 新增本土无症状感染者4390例 新增境外输入性新冠肺炎确诊病例2例 无新增境外输入性无症状感染者",
    "5月3日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年5月3日，新增本土新冠肺炎确诊病例260例 新增本土无症状感染者4722例 新增境外输入性新冠肺炎确诊病例3例 新增境外输入性无症状感染者1例",
    "5月2日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年5月2日，新增本土新冠肺炎确诊病例274例 新增本土无症状感染者5395例 无境外输入性新冠肺炎确诊病例 新增境外输入性无症状感染者1例",
    "5月1日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年5月1日，新增本土新冠肺炎确诊病例727例 新增本土无症状感染者6606例 新增境外输入性新冠肺炎确诊病例2例 新增境外输入性无症状感染者2例",
    "4月30日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年4月30日，新增本土新冠肺炎确诊病例788例 新增本土无症状感染者7084例 新增境外输入性新冠肺炎确诊病例1例 新增境外输入性无症状感染者2例",
    "4月29日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年4月29日，新增本土新冠肺炎确诊病例1249例 新增本土无症状感染者8932例 新增境外输入性新冠肺炎确诊病例1例 无新增境外输入性无症状感染者",
    "4月28日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年4月28日，新增本土新冠肺炎确诊病例5487例 新增本土无症状感染者9545例 新增境外输入性新冠肺炎确诊病例2例 新增境外输入性无症状感染者2例",
    "4月27日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "4月27日（0-24时）上海新增1292例本土新冠肺炎确诊病例，新增9330例本土无症状感染者",
    "4月26日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年4月26日，新增本土新冠肺炎确诊病例1606例 新增本土无症状感染者11956例 无新增境外输入性新冠肺炎确诊病例 无新增境外输入性无症状感染者",
    "【最新】4月25日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年4月25日，新增本土新冠肺炎确诊病例1661例 新增本土无症状感染者15319例 无新增境外输入性新冠肺炎确诊病例 无新增境外输入性无症状感染者",
    "4月24日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年4月24日，新增本土新冠肺炎确诊病例2472例 新增本土无症状感染者16983例 无新增境外输入性新冠肺炎确诊病例 新增境外输入性无症状感染者1例",
    "4月23日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年4月23日，新增本土新冠肺炎确诊病例1401例 新增本土无症状感染者19657例 无新增境外输入性新冠肺炎确诊病例 无新增境外输入性无症状感染者",
    "4月22日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年4月22日，新增本土新冠肺炎确诊病例2736例 新增本土无症状感染者20634例 无新增境外输入性新冠肺炎确诊病例 新增境外输入性无症状感染者1例",
    "4月21日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年4月21日，新增本土新冠肺炎确诊病例1931例 新增本土无症状感染者15698例 无新增境外输入性新冠肺炎确诊病例 无新增境外输入性无症状感染者",
    "4月20日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年4月20日，新增本土新冠肺炎确诊病例2634例 新增本土无症状感染者15861例 新增境外输入性新冠肺炎确诊病例1例 新增境外输入性无症状感染者2例",
    "4月19日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年4月19日，新增本土新冠肺炎确诊病例2494例 新增本土无症状感染者16407例 新增境外输入性新冠肺炎确诊病例1例 新增境外输入性无症状感染者1例",
    "4月18日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年4月18日，新增本土新冠肺炎确诊病例3084例 新增本土无症状感染者16358例 无新增境外输入性新冠肺炎确诊病例 新增境外输入性无症状感染者2例",  # 新增本土无症状感染者17332例 含既往无症状感染者转为确诊病例974例
    "4月17日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年4月17日，新增本土新冠肺炎确诊病例2417例 新增本土无症状感染者19831例 新增境外输入性新冠肺炎确诊病例3例 无新增境外输入性无症状感染者",
    "4月16日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年4月16日，新增本土新冠肺炎确诊病例3238例 新增本土无症状感染者21582例 新增境外输入性新冠肺炎确诊病例2例 无新增境外输入性无症状感染者",
    "4月15日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年4月15日，新增本土新冠肺炎确诊病例3590例 新增本土无症状感染者19923例 新增境外输入性新冠肺炎确诊病例4例 新增境外输入性无症状感染者2例",
    "4月14日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年4月14日，新增本土新冠肺炎确诊病例3200例 新增本土无症状感染者19872例 无新增境外输入性新冠肺炎确诊病例 新增境外输入性无症状感染者1例 治愈出院549例  解除医学观察8071例",
    "4月13日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年4月13日，新增本土新冠肺炎确诊病例2573例 新增本土无症状感染者25146例 无新增境外输入性新冠肺炎确诊病例 新增境外输入性无症状感染者1例 治愈出院741例  解除医学观察15406例",
    "【最新】4月12日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年4月12日，新增本土新冠肺炎确诊病例1189例 新增本土无症状感染者25141例 新增境外输入性新冠肺炎确诊病例1例 无新增境外输入性无症状感染者 治愈出院341例",
    "4月11日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年4月11日，新增本土新冠肺炎确诊病例994例 新增本土无症状感染者22348例 新增境外输入性新冠肺炎确诊病例4例 新增境外输入性无症状感染者1例 治愈出院391例",
    "4月10日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年4月10日，新增本土新冠肺炎确诊病例914例 新增本土无症状感染者25173例 新增境外输入性新冠肺炎确诊病例3例 治愈出院298例",
    "4月9日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年4月9日，新增本土新冠肺炎确诊病例1006例 新增本土无症状感染者23937例 新增境外输入性新冠肺炎确诊病例9例 新增境外输入性无症状感染者2例 治愈出院227例",
    "4月8日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年4月8日，新增本土新冠肺炎确诊病例1015例 新增本土无症状感染者22609例 新增境外输入性新冠肺炎确诊病例2例 新增境外输入性无症状感染者2例 治愈出院119例",
    "4月7日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年4月7日，新增本土新冠肺炎确诊病例824例 新增本土无症状感染者20398例 新增境外输入性新冠肺炎确诊病例4例 新增境外输入性无症状感染者3例 治愈出院57例",
    "4月6日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年4月6日，新增本土新冠肺炎确诊病例322例 新增本土无症状感染者19660例 新增境外输入性新冠肺炎确诊病例7例 新增境外输入性无症状感染者1例 治愈出院47例",
    "4月5日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年4月5日，新增本土新冠肺炎确诊病例311例 新增本土无症状感染者16766例 新增境外输入性新冠肺炎确诊病例4例 新增境外输入性无症状感染者1例 治愈出院41例",
    "4月4日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年4月4日，新增本土新冠肺炎确诊病例268例 新增本土无症状感染者13086例 新增境外输入性新冠肺炎确诊病例3例 新增境外输入性无症状感染者2例 治愈出院27例",
    "4月3日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年4月3日，新增本土新冠肺炎确诊病例425例 新增本土无症状感染者8581例 新增境外输入性新冠肺炎确诊病例8例 新增境外输入性无症状感染者4例 治愈出院30例",
    "4月2日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年4月2日，新增本土新冠肺炎确诊病例438例 新增本土无症状感染者7788例 新增境外输入性新冠肺炎确诊病例6例 新增境外输入性无症状感染者1例 治愈出院31例",
    "4月1日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年4月1日，新增本土新冠肺炎确诊病例260例 新增本土无症状感染者6051例 新增境外输入性新冠肺炎确诊病例2例 无新增境外输入性无症状感染者 治愈出院39例",
    "3月31日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年3月31日，新增本土新冠肺炎确诊病例358例 新增本土无症状感染者4144例 新增境外输入性新冠肺炎确诊病例7例 新增境外输入性无症状感染者1例 治愈出院26例",
    "3月30日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年3月30日，新增本土新冠肺炎确诊病例355例 新增本土无症状感染者5298例 新增境外输入性新冠肺炎确诊病例3例 新增境外输入性无症状感染者1例 治愈出院29例",
    "3月29日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年3月29日，新增本土新冠肺炎确诊病例326例 新增本土无症状感染者5656例 新增境外输入性新冠肺炎确诊病例3例 新增境外输入性无症状感染者2例 治愈出院45例",
    "3月28日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年3月28日，新增本土新冠肺炎确诊病例96例 新增本土无症状感染者4381例 新增境外输入性新冠肺炎确诊病例11例 新增境外输入性无症状感染者1例 治愈出院28例",
    "【最新】3月27日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年3月27日，新增本土新冠肺炎确诊病例50例 新增本土无症状感染者3450例 新增境外输入性新冠肺炎确诊病例10例 新增境外输入性无症状感染者4例 治愈出院37例",
    "3月26日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年3月26日，新增本土新冠肺炎确诊病例45例 新增本土无症状感染者2631例 新增境外输入性新冠肺炎确诊病例6例 新增境外输入性无症状感染者2例 治愈出院33例",
    "3月25日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年3月25日，新增本土新冠肺炎确诊病例38例 新增本土无症状感染者2231例 新增境外输入性新冠肺炎确诊病例9例 新增境外输入性无症状感染者2例 治愈出院21例  解除医学观察无症状感染者110例",
    "3月24日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年3月24日，新增本土新冠肺炎确诊病例29例 新增本土无症状感染者1580例 新增境外输入性新冠肺炎确诊病例12例 新增境外输入性无症状感染者5例 治愈出院63例  解除医学观察无症状感染者209例",
    "3月23日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年3月23日，新增本土新冠肺炎确诊病例4例 新增本土无症状感染者979例 新增境外输入性新冠肺炎确诊病例10例 新增境外输入性无症状感染者3例 治愈出院81例  解除医学观察无症状感染者230例",
    "3月22日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年3月22日，新增本土新冠肺炎确诊病例4例 新增本土无症状感染者977例 新增境外输入性新冠肺炎确诊病例10例 新增境外输入性无症状感染者4例 治愈出院23例  解除医学观察无症状感染者79例",
    "3月21日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年3月21日，新增本土新冠肺炎确诊病例31例 新增本土无症状感染者865例 新增境外输入性新冠肺炎确诊病例8例 新增境外输入性无症状感染者3例 治愈出院140例  解除医学观察无症状感染者186例",
    "3月20日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年3月20日，新增本土新冠肺炎确诊病例24例 新增本土无症状感染者734例 新增境外输入性新冠肺炎确诊病例17例 新增境外输入性无症状感染者2例 治愈出院44例  解除医学观察无症状感染者53例",
    "3月19日（0-24时）本市各区确诊病例、无症状感染者居住地信息",
    "上海2022年3月19日，新增本土新冠肺炎确诊病例17例 新增本土无症状感染者492例 新增境外输入性新冠肺炎确诊病例5例 新增境外输入性无症状感染者2例 治愈出院121例  解除医学观察无症状感染者68例",
    "上海2022年3月18日，新增本土新冠肺炎确诊病例8例 新增本土无症状感染者366例 新增境外输入性新冠肺炎确诊病例12例 新增境外输入性无症状感染者6例 治愈出院26例  解除医学观察无症状感染者15例",
    "上海2022年3月17日，新增本土新冠肺炎确诊病例57例 新增本土无症状感染者203例 新增境外输入性新冠肺炎确诊病例10例 新增境外输入性无症状感染者2例 治愈出院54例  解除医学观察无症状感染者38例",
    "上海2022年3月16日，新增本土新冠肺炎确诊病例8例 新增本土无症状感染者150例 新增境外输入性新冠肺炎确诊病例15例 新增境外输入性无症状感染者6例 治愈出院33例  解除医学观察无症状感染者19例",
    "“随申办”上线“居家健康监测证明”服务，部分区已可使用",
    "上海2022年3月15日，新增本土新冠肺炎确诊病例5例 新增本土无症状感染者197例 新增境外输入性新冠肺炎确诊病例10例 新增境外输入性无症状感染者1例 治愈出院33例  解除医学观察无症状感染者13例",
    "上海2022年3月14日，新增本土新冠肺炎确诊病例9例 新增本土无症状感染者130例 新增境外输入性新冠肺炎确诊病例12例 新增境外输入性无症状感染者3例 治愈出院15例  解除医学观察无症状感染者9例",
    "上海2022年3月13日，新增本土新冠肺炎确诊病例41例 新增本土无症状感染者128例 新增境外输入性新冠肺炎确诊病例16例 新增境外输入性无症状感染者2例 治愈出院33例  解除医学观察无症状感染者5例",
    "上海新增6例本土确诊病例和55例本土无症状感染者，两地列为中风险地区",
    "上海2022年3月12日，新增本土新冠肺炎确诊病例1例 新增本土无症状感染者64例 新增境外输入性新冠肺炎确诊病例13例 新增境外输入性无症状感染者1例 治愈出院26例  解除医学观察无症状感染者11例",
    "上海2022年3月11日，新增本土新冠肺炎确诊病例5例 新增本土无症状感染者78例 新增境外输入性新冠肺炎确诊病例22例 新增境外输入性无症状感染者9例 治愈出院20例  解除医学观察无症状感染者3例",
    "上海2022年3月10日，新增本土新冠肺炎确诊病例11例 新增本土无症状感染者64例 新增境外输入性新冠肺炎确诊病例32例 新增境外输入性无症状感染者10例 治愈出院17例  解除医学观察无症状感染者7例",
    "上海2022年3月9日，新增本土新冠肺炎确诊病例4例 新增本土无症状感染者76例 新增境外输入性新冠肺炎确诊病例42例 新增境外输入性无症状感染者16例 治愈出院18例  解除医学观察无症状感染者13例",
    "上海2022年3月8日，新增本土新冠肺炎确诊病例3例 新增本土无症状感染者62例 新增境外输入性新冠肺炎确诊病例26例 新增境外输入性无症状感染者10例 治愈出院32例  解除医学观察无症状感染者7例",
    "上海2022年3月7日，新增本土新冠肺炎确诊病例4例 新增本土无症状感染者51例 新增境外输入性新冠肺炎确诊病例36例 新增境外输入性无症状感染者10例 治愈出院27例  解除医学观察无症状感染者10例",
    "上海2022年3月6日，新增本土新冠肺炎确诊病例3例 新增本土无症状感染者45例 新增境外输入性新冠肺炎确诊病例32例 新增境外输入性无症状感染者16例 解除医学观察无症状感染者16例 治愈出院65例",
    "上海2022年3月5日，无新增本土新冠肺炎确诊病例 新增本土无症状感染者28例 新增境外输入性新冠肺炎确诊病例25例 新增境外输入性无症状感染者10例 解除医学观察无症状感染者4例 治愈出院8例",
    "上海2022年3月4日，新增本土新冠肺炎确诊病例3例 新增本土无症状感染者16例 新增境外输入性新冠肺炎确诊病例24例 新增境外输入性无症状感染者10例 治愈出院12例",
    "上海2022年3月3日，新增本土新冠肺炎确诊病例2例 新增本土无症状感染者14例 新增境外输入性新冠肺炎确诊病例43例 新增境外输入性无症状感染者21例 解除医学观察无症状感染者1例 治愈出院9例",
    "上海2022年3月2日，新增本土新冠肺炎确诊病例3例 新增本土无症状感染者5例 新增境外输入性新冠肺炎确诊病例39例 新增境外输入性无症状感染者19例 解除医学观察无症状感染者2例 治愈出院8例",
    "上海2022年3月1日，新增本土新冠肺炎确诊病例1例 新增本土无症状感染者1例 新增境外输入性新冠肺炎确诊病例37例 新增境外输入性无症状感染者17例 解除医学观察无症状感染者1例 治愈出院8例",
    "上海2022年2月28日，无新增本土新冠肺炎确诊病例 新增本土无症状感染者3例 新增境外输入性新冠肺炎确诊病例40例 新增境外输入性无症状感染者12例 治愈出院12例",
    "上海2022年2月27日，无新增本土新冠肺炎确诊病例 新增本土无症状感染者1例 新增境外输入性新冠肺炎确诊病例45例 新增境外输入性无症状感染者10例 治愈出院10例",
    "上海2022年2月26日，无新增本土新冠肺炎确诊病例 新增本土无症状感染者1例 新增境外输入性新冠肺炎确诊病例43例 新增境外输入性无症状感染者11例 治愈出院11例",
    "上海2022年2月25日，无新增本土新冠肺炎确诊病例 新增本土无症状感染者1例 新增境外输入性新冠肺炎确诊病例55例 新增境外输入性无症状感染者19例 治愈出院10例",
    "上海2022年2月24日，无新增本土新冠肺炎确诊病例 新增本土无症状感染者1例 新增境外输入性新冠肺炎确诊病例59例 新增境外输入性无症状感染者22例 治愈出院8例",
    "上海2022年2月23日，无新增本土新冠肺炎确诊病例 新增境外输入性新冠肺炎确诊病例24例 新增境外输入性无症状感染者8例 治愈出院8例",
    "上海2022年2月22日，无新增本土新冠肺炎确诊病例 新增境外输入性新冠肺炎确诊病例49例 新增境外输入性无症状感染者8例 治愈出院10例",
    "上海2022年2月21日，无新增本土新冠肺炎确诊病例 新增境外输入性新冠肺炎确诊病例32例 新增境外输入性无症状感染者1例 治愈出院7例",
    "上海2022年2月20日，无新增本土新冠肺炎确诊病例 新增本土无症状感染者2例 新增境外输入性新冠肺炎确诊病例27例 新增境外输入性无症状感染者3例 治愈出院3例",
    "上海2022年2月19日，无新增本土新冠肺炎确诊病例 新增境外输入23例  治愈出院13例",
    "上海2022年2月18日，无新增本土新冠肺炎确诊病例 新增境外输入性新冠肺炎确诊病例18例 新增境外输入性无症状感染者3例 治愈出院11例",
    "上海2022年2月17日，无新增本土新冠肺炎确诊病例 新增境外输入13例  治愈出院16例",
    "上海2022年2月16日，无新增本土新冠肺炎确诊病例 新增境外输入性新冠肺炎确诊病例9例 新增境外输入性无症状感染者1例 治愈出院19例",
]


def parse_line(line):
    if "无新增" in line:
        num = 0
    else:
        try:
            num = int(re.findall(r"\d+", line)[0])
        except Exception as e:
            print(e)
            print(line)
            num = 0
    if "解除" in line:
        return num, "relieve"
    elif "治愈" in line:
        return num, "cure"
    elif "境外" in line:
        if "无症状" in line:
            return num, "outbound_asymp"
        else:
            return num, "outbound_confirmed"
    else:
        if "无症状" in line:
            return num, "inbound_asymp"
        else:
            return num, "inbound_confirmed"


def parse_title(title_str):
    title_str = title_str.replace("  ", " ")
    if title_str.startswith("上海2022"):
        dic = {}
        date = title_str.split("，")[0][2:]
        raw = title_str.split("，")[1].split(" ")
        for r in raw:
            num_, type_ = parse_line(r)
            dic[type_] = num_
        return True, date, dic
    else:
        return False, None, None


df_dict = []
for l in p:
    flag, date, dic = parse_title(l)
    if flag:
        dic["date"] = date
        df_dict.append(dic)

df = pd.DataFrame(df_dict)
df["date"] = df["date"].apply(lambda x: dtm.datetime.strptime(x, "%Y年%m月%d日"))
df.fillna(0, inplace=True)
df.drop_duplicates(inplace=True)
df.sort_values(by=["date"], inplace=True)
df.to_csv("data/data_0515_SH.csv", index=False)
