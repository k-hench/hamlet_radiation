---
output: html_document
editor_options:
  chunk_output_type: console
css: highlight.css
---






# (git 2) Genotyping II (all callable sites)

This pipeline can be executed as follows:

```sh
cd $BASE_DIR/nf/02_genotyping_all_basepairs
source ../sh/nextflow_alias.sh
nf_run_allbp
```

## Summary

The genotyping procedure is controlled by the [**nextflow**](https://www.nextflow.io/) script `genotyping_all_basepairs.nf` (located under `$BASE_DIR/nf/02_genotyping_all_basepairs/`).
Based on an intermediate step from `genotyping.nf`, this script produces a data set that includes _all callable sites_  - that is SNPs as well a invariant sites that are covered by sequence.
Below is an overview of the steps involved in the genotyping process.
(The <span style="color:#4DAF4A">green dot</span> indicates the data input, <span style="color:#E41A1C">red arrows</span> depict output that is exported for further use.)

<div style="max-width:800px; margin:auto;">

```{=html}
<div id="htmlwidget-6cb2e5e2816767437751" style="width:672px;height:480px;" class="girafe html-widget"></div>
<script type="application/json" data-for="htmlwidget-6cb2e5e2816767437751">{"x":{"html":"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" id=\"svg_98855f50-f7c1-454a-8e25-73e30c18d8d5\" viewBox=\"0 0 648.00 648.00\">\n  <g>\n    <defs>\n      <clipPath id=\"svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_cl_1\">\n        <rect x=\"0.00\" y=\"0.00\" width=\"648.00\" height=\"648.00\"/>\n      <\/clipPath>\n    <\/defs>\n    <rect x=\"0.00\" y=\"0.00\" width=\"648.00\" height=\"648.00\" id=\"svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_el_1\" clip-path=\"url(#svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_cl_1)\" fill=\"#FFFFFF\" fill-opacity=\"1\" stroke-width=\"0.75\" stroke=\"#FFFFFF\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"round\"/>\n    <defs>\n      <clipPath id=\"svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_cl_2\">\n        <rect x=\"0.00\" y=\"0.00\" width=\"648.00\" height=\"648.00\"/>\n      <\/clipPath>\n    <\/defs>\n    <g clip-path=\"url(#svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_cl_2)\">\n      <text x=\"219.63\" y=\"433.35\" id=\"svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_el_2\" font-size=\"225.00pt\" font-weight=\"bold\" fill=\"#E0E0E0\" fill-opacity=\"1\" font-family=\"DejaVu Sans\">2<\/text>\n    <\/g>\n    <polyline points=\"542.44,37.35 542.11,39.70 541.75,42.27 541.38,44.83 541.02,47.39 540.66,49.95 540.29,52.51 539.93,55.08 539.57,57.64 539.20,60.20 538.84,62.76 538.48,65.32 538.11,67.89 537.75,70.45 537.39,73.01 537.03,75.57 536.66,78.14 536.30,80.70 535.94,83.26 535.57,85.82 535.21,88.38 534.85,90.95 534.48,93.51 534.12,96.07 533.76,98.63 533.39,101.20 533.03,103.76 532.67,106.32 532.31,108.88 531.94,111.44 531.58,114.01 531.22,116.57 530.85,119.13 530.49,121.69 530.13,124.25 529.76,126.82 529.40,129.38 529.04,131.94 528.67,134.50 528.31,137.07 527.95,139.63 527.59,142.19 527.22,144.75 526.86,147.31 526.50,149.88 526.13,152.44 525.77,155.00 525.41,157.56 525.04,160.13 524.68,162.69 524.32,165.25 523.95,167.81 523.59,170.37 523.23,172.94 522.86,175.50 522.50,178.06 522.14,180.62 521.78,183.18 521.41,185.75 521.05,188.31 520.69,190.87 520.32,193.43 519.96,196.00 519.60,198.56 519.23,201.12 518.87,203.68 518.51,206.24 518.14,208.81 517.78,211.37 517.42,213.93 517.06,216.49 516.69,219.06 516.33,221.62 515.97,224.18 515.60,226.74 515.24,229.30 514.88,231.87 514.51,234.43 514.15,236.99 513.79,239.55 513.42,242.11 513.06,244.68 512.70,247.24 512.34,249.80 511.97,252.36 511.61,254.93 511.25,257.49 510.88,260.05 510.52,262.61 510.16,265.17 509.79,267.74 509.43,270.30 509.07,272.86 508.73,275.22\" id=\"svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_el_3\" clip-path=\"url(#svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_cl_2)\" fill=\"none\" stroke-width=\"1.06698\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"butt\"/>\n    <polygon points=\"507.99,273.37 508.73,275.22 509.96,273.65\" id=\"svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_el_4\" clip-path=\"url(#svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_cl_2)\" fill=\"#000000\" fill-opacity=\"1\" stroke-width=\"1.06698\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"butt\"/>\n    <polyline points=\"613.82,427.31 612.94,426.12 611.82,424.60 610.70,423.08 609.58,421.56 608.46,420.04 607.34,418.51 606.22,416.99 605.10,415.47 603.98,413.95 602.86,412.43 601.74,410.91 600.62,409.39 599.50,407.87 598.38,406.34 597.26,404.82 596.14,403.30 595.01,401.78 593.89,400.26 592.77,398.74 591.65,397.22 590.53,395.69 589.41,394.17 588.29,392.65 587.17,391.13 586.05,389.61 584.93,388.09 583.81,386.57 582.69,385.04 581.57,383.52 580.45,382.00 579.33,380.48 578.21,378.96 577.09,377.44 575.97,375.92 574.85,374.39 573.73,372.87 572.60,371.35 571.48,369.83 570.36,368.31 569.24,366.79 568.12,365.27 567.00,363.74 565.88,362.22 564.76,360.70 563.64,359.18 562.52,357.66 561.40,356.14 560.28,354.62 559.16,353.09 558.04,351.57 556.92,350.05 555.80,348.53 554.68,347.01 553.56,345.49 552.44,343.97 551.32,342.44 550.19,340.92 549.07,339.40 547.95,337.88 546.83,336.36 545.71,334.84 544.59,333.32 543.47,331.79 542.35,330.27 541.23,328.75 540.11,327.23 538.99,325.71 537.87,324.19 536.75,322.67 535.63,321.14 534.51,319.62 533.39,318.10 532.27,316.58 531.15,315.06 530.03,313.54 528.90,312.02 527.78,310.49 526.66,308.97 525.54,307.45 524.42,305.93 523.30,304.41 522.18,302.89 521.06,301.37 519.94,299.84 518.82,298.32 517.70,296.80 516.58,295.28 515.46,293.76 514.34,292.24 513.22,290.72 512.34,289.53\" id=\"svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_el_5\" clip-path=\"url(#svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_cl_2)\" fill=\"none\" stroke-width=\"1.06698\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"butt\"/>\n    <polygon points=\"514.17,290.33 512.34,289.53 512.56,291.51\" id=\"svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_el_6\" clip-path=\"url(#svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_cl_2)\" fill=\"#000000\" fill-opacity=\"1\" stroke-width=\"1.06698\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"butt\"/>\n    <polyline points=\"501.09,287.69 500.07,288.40 498.81,289.28 497.56,290.17 496.30,291.05 495.04,291.93 493.78,292.81 492.53,293.70 491.27,294.58 490.01,295.46 488.75,296.34 487.50,297.22 486.24,298.11 484.98,298.99 483.72,299.87 482.47,300.75 481.21,301.63 479.95,302.52 478.69,303.40 477.43,304.28 476.18,305.16 474.92,306.05 473.66,306.93 472.40,307.81 471.15,308.69 469.89,309.57 468.63,310.46 467.37,311.34 466.12,312.22 464.86,313.10 463.60,313.99 462.34,314.87 461.09,315.75 459.83,316.63 458.57,317.51 457.31,318.40 456.06,319.28 454.80,320.16 453.54,321.04 452.28,321.92 451.03,322.81 449.77,323.69 448.51,324.57 447.25,325.45 446.00,326.34 444.74,327.22 443.48,328.10 442.22,328.98 440.97,329.86 439.71,330.75 438.45,331.63 437.19,332.51 435.94,333.39 434.68,334.27 433.42,335.16 432.16,336.04 430.91,336.92 429.65,337.80 428.39,338.69 427.13,339.57 425.88,340.45 424.62,341.33 423.36,342.21 422.10,343.10 420.85,343.98 419.59,344.86 418.33,345.74 417.07,346.63 415.82,347.51 414.56,348.39 413.30,349.27 412.04,350.15 410.79,351.04 409.53,351.92 408.27,352.80 407.01,353.68 405.76,354.56 404.50,355.45 403.24,356.33 401.98,357.21 400.73,358.09 399.47,358.98 398.21,359.86 396.95,360.74 395.70,361.62 394.44,362.50 393.18,363.39 391.92,364.27 390.67,365.15 389.65,365.87\" id=\"svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_el_7\" clip-path=\"url(#svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_cl_2)\" fill=\"none\" stroke-width=\"1.06698\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"butt\"/>\n    <polygon points=\"390.49,364.06 389.65,365.87 391.63,365.69\" id=\"svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_el_8\" clip-path=\"url(#svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_cl_2)\" fill=\"#000000\" fill-opacity=\"1\" stroke-width=\"1.06698\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"butt\"/>\n    <polyline points=\"376.60,375.02 375.91,375.50 374.71,376.34 373.51,377.19 372.31,378.03 371.11,378.87 369.90,379.72 368.70,380.56 367.50,381.40 366.30,382.24 365.10,383.09 363.90,383.93 362.70,384.77 361.49,385.62 360.29,386.46 359.09,387.30 357.89,388.14 356.69,388.99 355.49,389.83 354.28,390.67 353.08,391.52 351.88,392.36 350.68,393.20 349.48,394.04 348.28,394.89 347.08,395.73 345.87,396.57 344.67,397.42 343.47,398.26 342.27,399.10 341.07,399.94 339.87,400.79 338.66,401.63 337.46,402.47 336.26,403.32 335.06,404.16 333.86,405.00 332.66,405.84 331.46,406.69 330.25,407.53 329.05,408.37 327.85,409.22 326.65,410.06 325.45,410.90 324.25,411.75 323.04,412.59 321.84,413.43 320.64,414.27 319.44,415.12 318.24,415.96 317.04,416.80 315.84,417.65 314.63,418.49 313.43,419.33 312.23,420.17 311.03,421.02 309.83,421.86 308.63,422.70 307.42,423.55 306.22,424.39 305.02,425.23 303.82,426.07 302.62,426.92 301.42,427.76 300.22,428.60 299.01,429.45 297.81,430.29 296.61,431.13 295.41,431.97 294.21,432.82 293.01,433.66 291.80,434.50 290.60,435.35 289.40,436.19 288.20,437.03 287.00,437.87 285.80,438.72 284.60,439.56 283.39,440.40 282.19,441.25 280.99,442.09 279.79,442.93 278.59,443.78 277.39,444.62 276.18,445.46 274.98,446.30 273.78,447.15 272.58,447.99 271.38,448.83 270.69,449.31\" id=\"svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_el_9\" clip-path=\"url(#svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_cl_2)\" fill=\"none\" stroke-width=\"1.06698\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"butt\"/>\n    <polygon points=\"271.53,447.51 270.69,449.31 272.68,449.14\" id=\"svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_el_10\" clip-path=\"url(#svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_cl_2)\" fill=\"#000000\" fill-opacity=\"1\" stroke-width=\"1.06698\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"butt\"/>\n    <polyline points=\"257.64,458.47 257.04,458.89 255.85,459.72 254.67,460.56 253.48,461.39 252.29,462.22 251.10,463.06 249.91,463.89 248.73,464.72 247.54,465.56 246.35,466.39 245.16,467.22 243.98,468.06 242.79,468.89 241.60,469.72 240.41,470.56 239.22,471.39 238.04,472.22 236.85,473.06 235.66,473.89 234.47,474.72 233.28,475.56 232.10,476.39 230.91,477.22 229.72,478.06 228.53,478.89 227.34,479.72 226.16,480.56 224.97,481.39 223.78,482.22 222.59,483.06 221.41,483.89 220.22,484.72 219.03,485.56 217.84,486.39 216.65,487.22 215.47,488.06 214.28,488.89 213.09,489.72 211.90,490.56 210.71,491.39 209.53,492.22 208.34,493.06 207.15,493.89 205.96,494.72 204.77,495.56 203.59,496.39 202.40,497.22 201.21,498.06 200.02,498.89 198.84,499.72 197.65,500.56 196.46,501.39 195.27,502.22 194.08,503.06 192.90,503.89 191.71,504.72 190.52,505.56 189.33,506.39 188.14,507.22 186.96,508.06 185.77,508.89 184.58,509.72 183.39,510.56 182.20,511.39 181.02,512.22 179.83,513.06 178.64,513.89 177.45,514.72 176.27,515.56 175.08,516.39 173.89,517.22 172.70,518.06 171.51,518.89 170.33,519.72 169.14,520.56 167.95,521.39 166.76,522.22 165.57,523.06 164.39,523.89 163.20,524.72 162.01,525.56 160.82,526.39 159.63,527.22 158.45,528.06 157.26,528.89 156.07,529.72 154.88,530.56 153.70,531.39 153.09,531.81\" id=\"svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_el_11\" clip-path=\"url(#svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_cl_2)\" fill=\"none\" stroke-width=\"1.06698\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"butt\"/>\n    <polygon points=\"153.93,530.01 153.09,531.81 155.08,531.64\" id=\"svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_el_12\" clip-path=\"url(#svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_cl_2)\" fill=\"#000000\" fill-opacity=\"1\" stroke-width=\"1.06698\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"butt\"/>\n    <polyline points=\"140.04,540.97 139.47,541.37 138.29,542.20 137.10,543.03 135.92,543.86 134.74,544.69 133.56,545.52 132.37,546.35 131.19,547.18 130.01,548.01 128.82,548.84 127.64,549.67 126.46,550.50 125.27,551.33 124.09,552.16 122.91,552.99 121.73,553.82 120.54,554.65 119.36,555.48 118.18,556.31 116.99,557.14 115.81,557.97 114.63,558.80 113.44,559.63 112.26,560.46 111.08,561.28 109.90,562.11 108.71,562.94 107.53,563.77 106.35,564.60 105.16,565.43 103.98,566.26 102.80,567.09 101.62,567.92 100.43,568.75 99.25,569.58 98.07,570.41 96.88,571.24 95.70,572.07 94.52,572.90 93.33,573.73 92.15,574.56 90.97,575.39 89.79,576.22 88.60,577.05 87.42,577.88 86.24,578.71 85.05,579.54 83.87,580.37 82.69,581.20 81.50,582.03 80.32,582.86 79.14,583.69 77.96,584.52 76.77,585.35 75.59,586.18 74.41,587.01 73.22,587.84 72.04,588.67 70.86,589.50 69.68,590.33 68.49,591.16 67.31,591.99 66.13,592.82 64.94,593.65 63.76,594.48 62.58,595.31 61.39,596.14 60.21,596.97 59.03,597.80 57.85,598.63 56.66,599.46 55.48,600.29 54.30,601.12 53.11,601.95 51.93,602.78 50.75,603.61 49.56,604.44 48.38,605.27 47.20,606.10 46.02,606.93 44.83,607.76 43.65,608.59 42.47,609.42 41.28,610.25 40.10,611.08 38.92,611.91 37.74,612.74 36.55,613.57 35.98,613.97\" id=\"svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_el_13\" clip-path=\"url(#svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_cl_2)\" fill=\"none\" stroke-width=\"1.06698\" stroke=\"#E41A1C\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"butt\"/>\n    <polygon points=\"36.82,612.16 35.98,613.97 37.96,613.79\" id=\"svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_el_14\" clip-path=\"url(#svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_cl_2)\" fill=\"#E41A1C\" fill-opacity=\"1\" stroke-width=\"1.06698\" stroke=\"#E41A1C\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"butt\"/>\n    <g clip-path=\"url(#svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_cl_2)\">\n      <text transform=\"translate(532.38,187.56) rotate(-442)\" id=\"svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_el_15\" font-size=\"8.28pt\" font-family=\"DejaVu Sans\">vcf_cohort<\/text>\n    <\/g>\n    <g clip-path=\"url(#svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_cl_2)\">\n      <text transform=\"translate(537.37,342.25) rotate(-306)\" id=\"svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_el_16\" font-size=\"8.28pt\" font-family=\"DejaVu Sans\">ch_LG_ids<\/text>\n    <\/g>\n    <g clip-path=\"url(#svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_cl_2)\">\n      <text transform=\"translate(420.44,357.83) rotate(-395)\" id=\"svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_el_17\" font-size=\"8.28pt\" font-family=\"DejaVu Sans\">vcf_lg_combo<\/text>\n    <\/g>\n    <g clip-path=\"url(#svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_cl_2)\">\n      <text transform=\"translate(288.16,450.63) rotate(-395)\" id=\"svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_el_18\" font-size=\"8.28pt\" font-family=\"DejaVu Sans\">all_bp_by_location<\/text>\n    <\/g>\n    <g clip-path=\"url(#svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_cl_2)\">\n      <text transform=\"translate(178.06,527.87) rotate(-395)\" id=\"svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_el_19\" font-size=\"8.28pt\" font-family=\"DejaVu Sans\">all_bp_merged<\/text>\n    <\/g>\n    <g clip-path=\"url(#svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_cl_2)\">\n      <text transform=\"translate(64.94,607.22) rotate(-395)\" id=\"svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_el_20\" font-size=\"8.28pt\" font-family=\"DejaVu Sans\">filtered_snps<\/text>\n    <\/g>\n    <circle cx=\"543.56\" cy=\"29.45\" r=\"3.47pt\" id=\"svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_el_21\" clip-path=\"url(#svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_cl_2)\" fill=\"none\" stroke-width=\"0.708661\" stroke=\"#4DAF4A\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"round\"/>\n    <circle cx=\"618.55\" cy=\"433.73\" r=\"3.47pt\" id=\"svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_el_22\" clip-path=\"url(#svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_cl_2)\" fill=\"none\" stroke-width=\"0.708661\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"round\"/>\n    <circle cx=\"507.62\" cy=\"283.11\" r=\"3.47pt\" id=\"svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_el_23\" clip-path=\"url(#svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_cl_2)\" fill=\"none\" stroke-width=\"0.708661\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"round\"/>\n    <circle cx=\"383.12\" cy=\"370.44\" r=\"3.47pt\" id=\"svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_el_24\" clip-path=\"url(#svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_cl_2)\" fill=\"none\" stroke-width=\"0.708661\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"round\"/>\n    <circle cx=\"264.17\" cy=\"453.89\" r=\"3.47pt\" id=\"svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_el_25\" clip-path=\"url(#svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_cl_2)\" fill=\"none\" stroke-width=\"0.708661\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"round\"/>\n    <circle cx=\"146.57\" cy=\"536.39\" r=\"3.47pt\" id=\"svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_el_26\" clip-path=\"url(#svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_cl_2)\" fill=\"none\" stroke-width=\"0.708661\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"round\"/>\n    <circle cx=\"29.45\" cy=\"618.55\" r=\"3.47pt\" id=\"svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_el_27\" clip-path=\"url(#svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_cl_2)\" fill=\"none\" stroke-width=\"0.708661\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"round\"/>\n    <circle cx=\"543.56\" cy=\"29.45\" r=\"1.87pt\" id=\"svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_el_28\" clip-path=\"url(#svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_cl_2)\" fill=\"#4DAF4A\" fill-opacity=\"1\" stroke-width=\"0.708661\" stroke=\"#4DAF4A\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"round\" title=\"Channel.fromFilePairs\"/>\n    <circle cx=\"618.55\" cy=\"433.73\" r=\"1.87pt\" id=\"svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_el_29\" clip-path=\"url(#svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_cl_2)\" fill=\"#000000\" fill-opacity=\"1\" stroke-width=\"0.708661\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"round\" title=\"Channel.from\"/>\n    <circle cx=\"507.62\" cy=\"283.11\" r=\"1.87pt\" id=\"svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_el_30\" clip-path=\"url(#svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_cl_2)\" fill=\"#000000\" fill-opacity=\"1\" stroke-width=\"0.708661\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"round\" title=\"combine\"/>\n    <circle cx=\"383.12\" cy=\"370.44\" r=\"1.87pt\" id=\"svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_el_31\" clip-path=\"url(#svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_cl_2)\" fill=\"#000000\" fill-opacity=\"1\" stroke-width=\"0.708661\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"round\" title=\"joint_genotype_snps\"/>\n    <circle cx=\"264.17\" cy=\"453.89\" r=\"1.87pt\" id=\"svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_el_32\" clip-path=\"url(#svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_cl_2)\" fill=\"#000000\" fill-opacity=\"1\" stroke-width=\"0.708661\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"round\" title=\"merge_genotypes\"/>\n    <circle cx=\"146.57\" cy=\"536.39\" r=\"1.87pt\" id=\"svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_el_33\" clip-path=\"url(#svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_cl_2)\" fill=\"#000000\" fill-opacity=\"1\" stroke-width=\"0.708661\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"round\" title=\"filterSNPs\"/>\n    <circle cx=\"29.45\" cy=\"618.55\" r=\"1.87pt\" id=\"svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_el_34\" clip-path=\"url(#svg_98855f50-f7c1-454a-8e25-73e30c18d8d5_cl_2)\" fill=\"#000000\" fill-opacity=\"1\" stroke-width=\"0.708661\" stroke=\"#000000\" stroke-opacity=\"1\" stroke-linejoin=\"round\" stroke-linecap=\"round\" title=\"\"/>\n  <\/g>\n<\/svg>\n","js":null,"uid":"svg_98855f50-f7c1-454a-8e25-73e30c18d8d5","ratio":1,"settings":{"tooltip":{"css":" .tooltip_SVGID_ { padding:5px;background:black;color:white;border-radius:2px 2px 2px 2px ; position:absolute;pointer-events:none;z-index:999;}\n","offx":10,"offy":0,"use_cursor_pos":true,"opacity":0.9,"usefill":false,"usestroke":false,"delay":{"over":200,"out":500}},"hover":{"css":" .hover_SVGID_ { fill:orange;stroke:gray; }\n"},"hoverkey":{"css":" .hover_key_SVGID_ { stroke:red; }\n"},"hovertheme":{"css":" .hover_theme_SVGID_ { fill:green; }\n"},"zoom":{"min":1,"max":1},"capture":{"css":" .selected_SVGID_ { fill:red;stroke:gray; }\n","type":"multiple","only_shiny":true,"selected":[]},"capturekey":{"css":" .selected_key_SVGID_ { stroke:gray; }\n","type":"single","only_shiny":true,"selected":[]},"capturetheme":{"css":" .selected_theme_SVGID_ { stroke:gray; }\n","type":"single","only_shiny":true,"selected":[]},"toolbar":{"position":"topright","saveaspng":true},"sizing":{"rescale":true,"width":1}}},"evals":[],"jsHooks":[]}</script>
```
</div>

## Details of `genotyping_all_basepairs.nf`

### Data preparation

The nextflow script starts with a small header and then imports the joint genotyping likelihoods for all samples produced by `genotyping_all_basepairs.nf`.

<div class="kclass">

<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode">#<span class="hl opt">!/</span>usr<span class="hl opt">/</span>bin<span class="hl opt">/</span>env nextflow
<span class="hl slc">// git 2.1</span>
<span class="hl slc">// open genotype likelyhoods</span>
<span class="hl kwa">Channel</span>
	.fromFilePairs<span class="hl opt">(</span><span class="hl str">&quot;../../1_genotyping/1_gvcfs/cohort.g.vcf.{gz,gz.tbi}&quot;</span><span class="hl opt">)</span>
	.set<span class="hl opt">{</span> vcf_cohort <span class="hl opt">}</span>
</code>
</pre>
</div>

The genotyping of the different linkage groups is going to happen in parallel, so we need to initialize a channel for the 24 LGs.


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 2.2</span>
<span class="hl slc">// initialize LG channel</span>
<span class="hl kwa">Channel</span>
	.from<span class="hl opt">( (</span><span class="hl str">&#39;01&#39;</span>..<span class="hl str">&#39;09&#39;</span><span class="hl opt">) + (</span><span class="hl str">&#39;10&#39;</span>..<span class="hl str">&#39;19&#39;</span><span class="hl opt">) + (</span><span class="hl str">&#39;20&#39;</span>..<span class="hl str">&#39;24&#39;</span><span class="hl opt">) )</span>
	.set<span class="hl opt">{</span> ch_LG_ids <span class="hl opt">}</span>
</code>
</pre>
</div>

The genotyping likelihoods are combined, effectively linking the data set to the 24 parallel LGs.


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 2.3</span>
<span class="hl slc">// combine genotypes and LGs</span>
ch_LG_ids.combine<span class="hl opt">(</span> vcf_cohort <span class="hl opt">)</span>.set<span class="hl opt">{</span> vcf_lg_combo <span class="hl opt">}</span>
</code>
</pre>
</div>

The samples are jointly genotyped, independently for each LG and including invariant sites.


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 2.4</span>
<span class="hl slc">// actual genotyping step (including invariant sites)</span>
<span class="hl kwa">process</span> joint_genotype_snps <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&quot;L_O88g90h_LGs_genotype&quot;</span>

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> lg <span class="hl opt">),</span> vcfId<span class="hl opt">,</span> <span class="hl kwc">file</span><span class="hl opt">(</span> vcf <span class="hl opt">)</span> <span class="hl kwa">from</span> vcf_lg_combo

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> <span class="hl str">&#39;all&#39;</span> <span class="hl opt">),</span> <span class="hl kwc">val</span><span class="hl opt">(</span> lg <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;all_site*.vcf.gz&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;all_site*.vcf.gz.tbi&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> all_bp_by_location

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	gatk --java-options &quot;-Xmx85g&quot; \</span>
<span class="hl str">		GenotypeGVCFs \</span>
<span class="hl str">		-R=\$BASE_DIR/ressources/HP_genome_unmasked_01.fa \</span>
<span class="hl str">		-L=LG</span><span class="hl ipl">${lg}</span> <span class="hl str">\</span>
<span class="hl str">		-V=</span><span class="hl ipl">${vcf[0]}</span> <span class="hl str">\</span>
<span class="hl str">		-O=intermediate.vcf.gz \</span>
<span class="hl str">		--include-non-variant-sites=true</span>
<span class="hl str"></span>
<span class="hl str">	gatk --java-options &quot;-Xmx85G&quot; \</span>
<span class="hl str">		SelectVariants \</span>
<span class="hl str">		-R=\$BASE_DIR/ressources/HP_genome_unmasked_01.fa \</span>
<span class="hl str">		-V=intermediate.vcf.gz \</span>
<span class="hl str">		--select-type-to-exclude=INDEL \</span>
<span class="hl str">		-O=all_sites.LG</span><span class="hl ipl">${lg}</span><span class="hl str">.vcf.gz</span>
<span class="hl str"></span>
<span class="hl str">	rm intermediate.*</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>

The genotypes of the different LGs are merged.


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 2.5</span>
<span class="hl slc">// merge all LGs</span>
<span class="hl kwa">process</span> merge_genotypes <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&#39;L_78g5h_merge_genotypes&#39;</span>
	<span class="hl kwb">echo</span> true

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">val</span><span class="hl opt">(</span> dummy <span class="hl opt">),</span>  <span class="hl kwc">val</span><span class="hl opt">(</span> lg <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> vcf <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> tbi <span class="hl opt">)</span> <span class="hl kwa">from</span> all_bp_by_location.groupTuple<span class="hl opt">()</span>

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;all_sites.vcf.gz&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> all_bp_merged

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	INPUT=\$(ls -1 *vcf.gz | sed &#39;s/^/ -I /g&#39; | cat \$( echo ))</span>
<span class="hl str"></span>
<span class="hl str">	gatk --java-options &quot;-Xmx85g&quot; \</span>
<span class="hl str">		GatherVcfs \</span>
<span class="hl str">		\$INPUT \</span>
<span class="hl str">		-O=all_sites.vcf.gz</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>

The genotypes are hard filtered based on various genotyping scores.


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 2.6</span>
<span class="hl slc">// quality based filtering</span>
<span class="hl kwa">process</span> filterSNP_first <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&#39;L_105g30h_filter_gt1&#39;</span>

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwc">file</span><span class="hl opt">(</span> vcf <span class="hl opt">)</span> <span class="hl kwa">from</span> all_bp_merged

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;intermediate.filterd.vcf.gz&quot;</span> <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;intermediate.filterd.vcf.gz.tbi&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> filtered_snps_first

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	module load openssl1.0.2</span>
<span class="hl str"></span>
<span class="hl str">	tabix -p vcf</span> <span class="hl ipl">${vcf}</span>
<span class="hl str"></span>
<span class="hl str">	gatk --java-options &quot;-Xmx75G&quot; \</span>
<span class="hl str">		VariantFiltration \</span>
<span class="hl str">		-R=\$BASE_DIR/ressources/HP_genome_unmasked_01.fa \</span>
<span class="hl str">		-V</span> <span class="hl ipl">${vcf}</span> <span class="hl str">\</span>
<span class="hl str">		-O=intermediate.vcf.gz \</span>
<span class="hl str">		--filter-expression &quot;QD &lt; 2.5&quot; \</span>
<span class="hl str">		--filter-name &quot;filter_QD&quot; \</span>
<span class="hl str">		--filter-expression &quot;FS &gt; 25.0&quot; \</span>
<span class="hl str">		--filter-name &quot;filter_FS&quot; \</span>
<span class="hl str">		--filter-expression &quot;MQ &lt; 52.0 || MQ &gt; 65.0&quot; \</span>
<span class="hl str">		--filter-name &quot;filter_MQ&quot; \</span>
<span class="hl str">		--filter-expression &quot;MQRankSum &lt; -0.2 || MQRankSum &gt; 0.2&quot; \</span>
<span class="hl str">		--filter-name &quot;filter_MQRankSum&quot; \</span>
<span class="hl str">		--filter-expression &quot;ReadPosRankSum &lt; -2.0 || ReadPosRankSum &gt; 2.0 &quot; \</span>
<span class="hl str">		--filter-name &quot;filter_ReadPosRankSum&quot; \</span>
<span class="hl str">		--QUIET true &amp;&gt; var_filt.log</span>
<span class="hl str"></span>
<span class="hl str">	gatk --java-options &quot;-Xmx75G&quot; \</span>
<span class="hl str">		SelectVariants \</span>
<span class="hl str">		-R=\$BASE_DIR/ressources/HP_genome_unmasked_01.fa \</span>
<span class="hl str">		-V=intermediate.vcf.gz \</span>
<span class="hl str">		-O=intermediate.filterd.vcf.gz \</span>
<span class="hl str">		--exclude-filtered \</span>
<span class="hl str">		--QUIET true \</span>
<span class="hl str">		--verbosity ERROR  &amp;&gt; var_select.log</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>

A second filtering is based on the missingness of samples.


<div class="sourceCode">
<pre class="sourceCode">
<code class="sourceCode"><span class="hl slc">// git 2.7</span>
<span class="hl slc">// missingness based filtering</span>
<span class="hl slc">// the resulting vcf file represents</span>
<span class="hl slc">// the &#39;all BP&#39; data set</span>
<span class="hl kwa">process</span> filterSNP_second <span class="hl opt">{</span>
	<span class="hl kwb">label</span> <span class="hl str">&#39;L_105g30h_filter_gt2&#39;</span>
	<span class="hl kwb">publishDir</span> <span class="hl str">&quot;../../1_genotyping/3_gatk_filtered/&quot;</span><span class="hl opt">,</span> mode<span class="hl opt">:</span> <span class="hl str">&#39;copy&#39;</span>

	<span class="hl kwb">input</span><span class="hl opt">:</span>
	<span class="hl kwa">set</span> <span class="hl kwc">file</span><span class="hl opt">(</span> vcf <span class="hl opt">),</span> <span class="hl kwc">file</span><span class="hl opt">(</span> tbi <span class="hl opt">)</span> <span class="hl kwa">from</span> filtered_snps_first

	<span class="hl kwb">output</span><span class="hl opt">:</span>
	<span class="hl kwc">file</span><span class="hl opt">(</span> <span class="hl str">&quot;filterd.allBP.vcf.gz&quot;</span> <span class="hl opt">)</span> <span class="hl kwa">into</span> filtered_snps

	<span class="hl kwb">script</span><span class="hl opt">:</span>
	<span class="hl str">&quot;&quot;&quot;</span>
<span class="hl str">	module load openssl1.0.2</span>
<span class="hl str"></span>
<span class="hl str">	vcftools \</span>
<span class="hl str">		--gzvcf</span> <span class="hl ipl">${vcf}</span> <span class="hl str">\</span>
<span class="hl str">		--max-missing-count 17 \</span>
<span class="hl str">		--stdout  \</span>
<span class="hl str">		--recode | \</span>
<span class="hl str">		bgzip &gt; filterd.allBP.vcf.gz</span>
<span class="hl str">	&quot;&quot;&quot;</span>
<span class="hl opt">}</span>
</code>
</pre>
</div>
</div>

Finally, we are done with the second version of genotyping.

---
