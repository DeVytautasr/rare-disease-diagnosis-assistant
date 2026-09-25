# Naudojimas

Šis dokumentas paaiškina, **ką rodo įrankis ir kaip tai skaityti**.
Diegimas aprašytas atskirai — žr. DIEGIMAS.md.

---

## Svarbiausia iš karto: ką šis įrankis daro ir ko nedaro

Įrankis **skaito jau paruoštą kandidatų rinkinį** (`.bcf` arba `.vcf` failą) ir
parodo, kokie duomenys BAM faile tuos kandidatus pagrindžia arba nepagrindžia.

**Jis pats kandidatų neieško.** Kandidatų rinkinį pagamina atskira programa —
`delly` — ir tai yra visiškai atskiras žingsnis, kurio čia nėra ir kuris
**nebuvo įdiegtas** kartu su šiuo įrankiu.

> Jeigu skaitote tik šį dokumentą, nesusidarykite įspūdžio, kad įdiegėte visą
> analizės grandinę. Įdiegėte **antrąją jos dalį** — skaitytuvą.

### Ką kainuoja pirmoji dalis (kandidatų gaminimas)

`delly` reikia dalykų, kurių šiam įrankiui nereikia: **referencinio genomo**
(apie 3 GB), pačios `delly` programos ir nemažai skaičiavimo laiko.

Pavyzdinė komanda (taip buvo pagamintas demonstracinis rinkinys):

```bash
delly sr \
  -g GRCh38_full_analysis_set_plus_decoy_hla.fa \
  -x human.hg38.excl.tsv \
  -o meginys.bcf \
  --threads 1 \
  meginys.bam
```

Kiek tai trunka — išmatuota:

| BAM dydis | Trukmė | Didžiausias atminties poreikis |
|---|---|---|
| ~1,6 GB (dvi chromosomos) | **39 s** | mažas |
| ~40 GB (visas genomas) | **46 min** ir **2 val. 4 min** | apie **1,2 GB** |

Tai daroma **vieną kartą** kiekvienam mėginiui. Gautą `.bcf` failą (paprastai
apie 130 KB) po to galima skaityti šiuo įrankiu kiek nori kartų.

---

## Darbo eiga

### 1. Įkelti kandidatų rinkinį

Viršuje pasirinkite kandidatų rinkinį iš sąrašo ir paspauskite **Load**.
Sąraše matomos **etiketės**, ne failų keliai.

### 2. Filtrų grandinė — ką pašalina kiekvienas žingsnis

Lentelėje kiekviena eilutė — vienas filtras. Stulpeliai:

| Stulpelis | Ką reiškia |
|---|---|
| **step** | filtro pavadinimas |
| **cut-off** | riba, nuo kurios jis veikia |
| **where the cut-off came from** | iš kur ta riba paimta |
| **still remaining** | kiek kandidatų liko **po** šio žingsnio |
| **removed by this step** | kiek pašalino **šioje vietoje** grandinėje |
| **this step would remove on its own** | kiek būtų pašalinęs **vienas**, be kitų |

Du paskutiniai stulpeliai dažnai skiriasi, ir tai svarbu. Jei filtras
„pašalino 0“, bet „vienas būtų pašalinęs 300“, jis **nėra nenaudingas** —
tiesiog ankstesni filtrai tuos įrašus jau buvo pašalinę. Įrankis tokius
žingsnius pažymi atskirai, kad nepasirodytų, jog filtras nieko nedaro.

Eilė turi reikšmės. Tie patys filtrai kita tvarka duoda tuos pačius galutinius
kandidatus, bet kitokius tarpinius skaičius.

**Jei žingsnis neatliktas**, virš lentelės atsiras raudonas užrašas su
paaiškinimu. Dažniausia priežastis — nerastas delly neįtraukiamų sričių
šablonas. Tokiu atveju grandinė veikia, tik be to vieno žingsnio, ir apie tai
pasako garsiai.

### 3. Atverti kandidatą

Paspauskite kandidato eilutę. Įrankis paleis keturis įrodymų sluoksnius abiejose
lūžio taško pusėse.

### 4. Įvesti koordinatę ranka

Laukeliuose įrašykite chromosomą ir poziciją ir paspauskite **Assess**.
Koordinatė **nebūtinai turi būti iš kandidatų rinkinio** — galima tikrinti bet
kurią vietą.

Atsakyme visada matysite `position_provenance` — iš kur ta pozicija atsirado.
Jei įvedėte ranka, ten bus parašyta, kad jos **nepatvirtino joks kandidatų
rinkinys**. Tai apsauga nuo savęs apgaudinėjimo: radus „įrodymų“ ranka įvestoje
vietoje, tai dar nereiškia, kad ten yra tikras lūžis.

### 5. Palyginti du rinkinius

Skydelyje **Compare** pasirinkite du kandidatų rinkinius ir toleranciją
bazėmis (pvz. 500). Įrankis parodys, kiek sutampa ir kiek yra tik viename.

Naudinga, kai tas pats mėginys apdorotas skirtingais nustatymais — matyti,
ką pakeitimas realiai pridėjo arba atėmė.

---

## Keturi įrodymų sluoksniai

| Sluoksnis | Ką skaičiuoja | Ką reiškia |
|---|---|---|
| **Poriniai skaitiniai** (discordant pairs) | poras, kurių du galai nusėdo ne ten, kur turėtų | rodo, su kuo sujungta |
| **Nukirsti skaitiniai** (soft-clipped) | skaitinius, kurių galas „nukerpamas“ toje pačioje vietoje | rodo **tikslią** lūžio vietą |
| **Perskelti skaitiniai** (split reads) | skaitinius, kurių dalis nusėdo kitur | pats tiesiausias įrodymas |
| **Skaitymo gylis** (read depth) | ar padaugėjo/sumažėjo medžiagos | rodo iškritas ir dublikacijas |

Kiekvienas vertinamas atskirai nuo 0 iki 25, iš viso 0–100.
Vertinimas: **70 ir daugiau — „strong“**, **40 ir daugiau — „moderate“**,
daugiau nei 0 — „weak“.

Prie kiekvieno skaičiaus rodomas mygtukas, atidarantis **tikslų įrankio
atsakymą**. Nė vienas ekrane matomas skaičius neatsirado kitaip.

---

## „QUALITY-LIMITED“ — tai ne blogas įvertinimas

Jei per daug skaitinių toje vietoje nusėdo nepatikimai (daugiau nei 40 %
žemiau MAPQ 20), įrankis **neskaičiuoja** bendro įvertinimo ir parašo
`QUALITY-LIMITED`.

> **Tai reiškia „nepamatuota“, o ne „pamatuota ir mažai“.**

Skirtumas esminis. „weak“ reiškia: pažiūrėjome ir beveik nieko neradome.
`QUALITY-LIMITED` reiškia: **negalėjome pažiūrėti** — ši genomo vieta tokia
pasikartojanti, kad nežinia, ar skaitiniai apskritai iš čia.

Keturi atskiri sluoksniai vis tiek rodomi. Skaitykite juos.

---

## Kodėl subalansuota translokacija niekada negaus „strong“

Tai aritmetika, ne duomenų trūkumas.

Subalansuotoje translokacijoje **medžiagos nei pridedama, nei atimama**.
Todėl skaitymo gylio sluoksnis teisingai duoda **0 taškų iš 25** — ir tai
teisingas atsakymas, ne klaida.

Antra: jei pertvarkyta **tik viena iš dviejų** chromosomos kopijų
(heterozigotinis atvejis), maždaug **pusė** skaitinių toje vietoje ateina iš
**sveikosios kopijos**. Todėl porinių skaitinių dalis niekada nepasiekia 0,5,
kurios reikalauja aukščiausia to sluoksnio pakopa.

Sudėjus: du sluoksniai iš keturių yra apriboti iš anksto, ir aukščiausia
juosta pasidaro **nepasiekiama**.

Įrankis tai apskaičiuoja **kiekvienai vietai atskirai** ir parodo laukuose
`attainable_here` (didžiausias čia pasiekiamas įvertinimas) ir `strong_band`
(nuo kiek prasideda „strong“). Jei pirmasis mažesnis už antrąjį — „strong“
čia neįmanomas, kad ir kokie geri būtų duomenys.

> Subalansuotą translokaciją vertinkite pagal **keturis atskirus matavimus**,
> o ne pagal juostą, į kurią pateko bendras skaičius.

---

## Demonstracinio rinkinio ribos

Demonstracinis rinkinys yra **1,6 MB**. Jame yra du mėginiai, o jų BAM failai
**iškarpyti**: palikti tik reikalingi ruožai, o ne visos chromosomos. Taip
daroma sąmoningai — pilni failai svertų po 1,6 GB kiekvienas.

| Mėginys | Kas jame palikta | Ką parodo |
|---|---|---|
| `DEMO_CLEAN` | `chr20:200 000` ir `chr21:14 100 000`, po ±10 000 bazių | trys iš keturių sluoksnių; įvertis 40,0/100 „moderate“; aukščiausia juosta čia nepasiekiama (daugiausia 57,5, reikia 70) |
| `DEMO_REPEAT` | `chr20:25 800 000` ir `chr21:7 600 000`, po ±10 000 bazių | `QUALITY-LIMITED` — įvertis sulaikytas |

> *Pataisyta 2026-09-25.* Čia buvo parašyta „visi keturi sluoksniai; įvertis
> 47,5/100“ — tai pirmojo implantų rinkinio reikšmė; to rinkinio įrašai prarasti.
> Rinkinys perkurtas, ir dabar `make_demo_bundle.py` šias reikšmes kiekvieną kartą
> išmatuoja pačiame iškarpytame faile ir įrašo į `DEMO.md`.

| Veikia visiškai | Veikia tik ruožuose |
|---|---|
| filtrų grandinė (896 sujungimai po sulydymo) | keturi įrodymų sluoksniai |
| dviejų rinkinių palyginimas | ranka įvesta koordinatė |

Filtrų grandinė ir palyginimas veikia **pilnai**, nes jie skaito tik `.bcf`
failus, o tie nesukarpyti.

Įvedus koordinatę už ruožų ribų, sluoksniai parodys, kad **duomenų ten nėra**
(`assessable: false`) — tai ne klaida, o sąžiningas atsakymas. Tikruose
duomenyse tokio apribojimo nėra.

> **Patikrinta:** tuose ruožuose iškarpytas failas duoda **tą patį** įvertį kaip
> ir pilnas — 40,0/100 (pirmajame implantų rinkinyje buvo 47,5/100). Tai nėra savaime suprantama: pirmasis bandymas davė
> 43,3, nes sluoksnių tinkamumas nustatomas pagal failo pradžios skaitinius, o
> iškarpytame faile jie yra kitokie. Tai ištaisyta.

---

## Ko šis įrankis pasakyti negali

**Jis tikrina pozicijas, bet jų neieško.** Kiekviena čia tikrinama vieta arba
atėjo iš kandidatų rinkinio, arba buvo įvesta ranka. Tikras lūžis, kurio
neaptiko `delly` ir kurio niekas neįvedė, čia **niekada nepasirodys**.

**Kontroliuotame teste rasti 16 iš 24 žinomų lūžio taškų.** Švarioje,
vienareikšmiškai skaitomoje sekoje — 8 iš 8. Šalia pasikartojančių sričių —
8 iš 8. Ten, kur seka skaitoma dviprasmiškai — **0 iš 8**. Visi praleisti buvo
prarasti kandidatų paieškos žingsnyje, ne filtruose.

> *Pataisyta 2026-09-25.* Čia buvo parašyta „14 iš 24“ ir „šalia pasikartojančių
> sričių 6 iš 8“ — tai pirmojo bandymo rezultatas; jo įrašai prarasti. Testas
> perkurtas kaip naujas eksperimentas. Skirtumą lemia IMP06: tada praleistas,
> dabar aptiktas; perlyginus jo skaitinius be `.alt` failo jis vėl prarandamas,
> taigi skirtumą paaiškina ALT-aware lyginimas, o ne atsitiktinė imtis.

**14 iš 16 ribų yra autoriaus sprendimas, o ne kalibruotos vertės.** Jos
nepatikrintos su patvirtintų teigiamų ir neigiamų atvejų rinkiniu. Prie
kiekvienos ribos parašyta, iš kur ji. „Autoriaus sprendimas“ skaitykite taip:
protingas žmogus būtų pasirinkęs kitaip, ir rezultatas būtų kitoks.

**Skaitymo gylio matavimas prie tokio padengimo nepatikimas.** Jo riba nustatyta
su maždaug dešimt kartų didesnio padengimo duomenimis. Vertinkite jį kaip
silpną papildomą požymį, ne daugiau.

**Įvertinimas nėra tikimybė.** Tai suskaidytas aprašas, kas suveikė, o ne
apskaičiuota tikimybė, kad variantas tikras.

**Tai magistro darbo prototipas, ne klinikinis įrankis.** Nė vienas čia matomas
skaičius netinka klinikiniam sprendimui be nepriklausomo patvirtinimo.
