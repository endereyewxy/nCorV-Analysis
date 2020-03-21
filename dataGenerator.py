#!/bin/python3

# 数据来源：国家卫生健康委员会官方网站 http://www.nhc.gov.cn/xcs/yqtb/list_gzbd.shtml

import requests, csv, re, lxml.etree

HOST = 'http://www.nhc.gov.cn'

# 由于官网的坑爹验证机制，每次运行脚本必须从浏览器手动更新cookies
COOKIES = """
oHAcoULcWCQb80S=DaJBaNEDHUOY9Sx53rEPCYPlsYxumABpXEyPPsFhjyvvNwS0vS_DkJC9N365tWN4; insert_cookie=96816998; oHAcoULcWCQb80T=4Ss9Xi.A_uGdK_GFS4MU0Iuz4F3VctSw_IDNB.ugAulSbNl6OtUm1OMsNINoE6VTPnG5xV39Pge7A0vXRJg2DMNTLftQ1ahA700RsQS_qpcYuU..ARGP7eesCL9Gp1WEy7vTY1Xx5uNR4s_Hh7yzeRC7L3a9AewHO.9OWn2fu8VFBbdoOkc4WGTcNeYJNrRfESSo5ggKN0MB8Jh6ziSg16VDFieKuOOBJaRXfALzj6RYx21w1TWJdQVY2qILfMjLyIkZoHNvpxPrlnmf_9.6oz6GDZUBRctpS7MzwHi95HftuhWsWuk8EJwYpZwgSGG6evaW; security_session_verify=90a15877f97758e7e6d37eae0d2ab7db; yfx_c_g_u_id_10006654=_ck20032112162410444335353995809; yfx_f_l_v_t_10006654=f_t_1584764184042__r_t_1584764184042__v_t_1584772110384__r_c_0
"""
HEADERS = {
    'Host': 'www.nhc.gov.cn',
    'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64; rv:74.0) Gecko/20100101 Firefox/74.0',
    'Cookie': COOKIES.strip()
}
# 辣鸡卫建委，天天改报告格式
PATTERNS = [r'新增确诊\D*(\d+)例', r'新增疑似\D*(\d+)例', r'新增治愈\D*(\d+)例', r'新增死亡\D*(\d+)例', r'[现|共]有确诊\D*(\d+)例', r'[现|共]有疑似\D*(\d+)例', '累计\D*确诊\D*(\d+)例', r'累计\D*治愈\D*(\d+)例', r'累计\D*死亡\D*(\d+)例']
TITLE = ['日期', '新增确诊', '新增疑似', '新增治愈', '新增死亡', '现有确诊', '现有疑似', '累计确诊', '累计治愈', '累计死亡']

def getHtml(url):
    resp = requests.get(url, headers=HEADERS)
    if resp.status_code != 200:
        print('网络连接错误：', resp.status_code)
        return None
    return lxml.etree.HTML(resp.text)

def main():
    url_indices = [HOST + '/xcs/yqtb/list_gzbd' + suffix + '.shtml' for suffix in ['', '_2', '_3']]
    url_targets = []
    for url in url_indices:
        root = getHtml(url)
        if root is None:
            return
        for link in root.xpath('//li'):
            date = link.xpath('span')[0].text
            # 由于格式不统一，忽略2020-02-07之前的数据
            if date[5:7] == '01' or (date[5:7] == '02' and int(date[8:]) < 7):
                continue
            url_targets.append((date, HOST + link.xpath('a')[0].get('href')))
    data = []
    for date, url in url_targets:
        page = getHtml(url)
        if page is None:
            return
        content = ''
        for item in page.xpath('//div[@id]')[0].xpath('p')[0].iter():
            content += str(item.text if item.tag == 'p' else item.tail)
        row = [date]
        for pat in PATTERNS:
            row.append(re.search(pat, content).group(1))
        print(row)
        data.append(row)
    with open('data.csv', 'w') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(TITLE)
        for row in data:
            writer.writerow(row)

if __name__ == '__main__':
    main()

