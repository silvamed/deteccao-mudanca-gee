// Deteccao de mudanca de cobertura em propriedades militares (Mata Atlantica)
// App para o Code Editor do Google Earth Engine - versao 2.7
// PFC de Engenharia Cartografica, IME - Pedro Medeiros e Pablo Artte
//
// Como usar:
//   1. Desenhe o poligono da area no mapa da esquerda ou cole o id de um asset
//      e clique em "Usar asset".
//   2. Clique em "Verificar bioma".
//   3. Clique em "Processar". A amostra de validacao (SHP de pontos e CSV com a
//      area dos estratos) vai direto para a aba Tasks.
//   4. Use os botoes de "Saidas" para exportar os produtos para o Drive.
//
// Os parametros so foram testados na Mata Atlantica. Fora dela o app pede
// confirmacao antes de processar.
// No app publicado as exportacoes e o carregamento por asset nao funcionam;
// rode pelo Code Editor.
// O limite da propriedade nao fica salvo no codigo: use asset privado ou desenhe.


// ---------------------------------------------------------------------------
// Configuracao
// ---------------------------------------------------------------------------
var ANO_INI = 2019;
var ANO_FIM = 2025;

var BANDAS    = ['B2','B3','B4','B8','B11','B12'];
var CS_LIMIAR = 0.55;    // Cloud Score+ (cs_cdf)
var ESCALA    = 10;      // m, exportacoes e tabelas
var ESCALA_UI = 30;      // m, numeros do painel (so para leitura rapida)
var MARGEM    = 2000;    // m em volta da area, para as composicoes exportadas
var PASTA_EXPORT = 'PFC_degradacao';
var VERSAO = '2.7';
var SEM_DADO = 255;      // valor dos rasters de classe fora da area ou sem dado

var N_VALIDACAO_POR_ESTRATO = 100;
var SEMENTE_AMOSTRA = 42;
var SEMENTE_ORDEM   = 7;

var S2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED');
var CS = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED');
var CHIRPS = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY');

// Janela jun-set, limiares e endmembers [B2, B3, B4, B8, B11, B12] em
// reflectancia. Endmembers = media de pixels puros tirados da propria imagem.
var PARAM = {
  rotulo: 'Mata Atlantica', sigla: 'MA',
  mes_ini: 6, mes_fim: 9,
  lim_solo: 0.30, lim_ndvi: 0.45,
  em_veg:    [0.033, 0.057, 0.040, 0.274, 0.182, 0.095],
  em_solo:   [0.087, 0.136, 0.186, 0.280, 0.334, 0.255],
  em_sombra: [0.010, 0.010, 0.010, 0.012, 0.010, 0.010],
  nota: 'Parametros testados apenas na Mata Atlantica.'
};

var ROTULO_BIOMA = {
  MATA_ATLANTICA: 'Mata Atlantica', CERRADO: 'Cerrado', CAATINGA: 'Caatinga',
  PAMPA: 'Pampa', PANTANAL: 'Pantanal', AMAZONIA: 'Amazonia'
};


// ---------------------------------------------------------------------------
// Bioma (RESOLVE Ecoregions 2017)
// ---------------------------------------------------------------------------
// Usa a ecorregiao com maior area dentro do poligono. A conversao para bioma
// e feita pelo nome, porque Amazonia e Mata Atlantica estao no mesmo bioma da
// WWF. Na duvida, confira no mapa de biomas do IBGE.
var ECO_MATA_ATLANTICA = [
  'serra do mar coastal forests',
  'bahia coastal forests',
  'bahia interior forests',
  'pernambuco coastal forests',
  'pernambuco interior forests',
  'alto parana atlantic forests',
  'araucaria moist forests',
  'atlantic coast restingas',
  'brazilian atlantic dry forests',
  'caatinga enclaves moist forests'
];

// Minusculas e sem acento ("Alto Parana" vem com acento no RESOLVE).
function normalizar(txt) {
  var s = String(txt || '').toLowerCase();
  if (s.normalize) { s = s.normalize('NFD').replace(/[\u0300-\u036f]/g, ''); }
  return s;
}

function biomaDaEcorregiao(ecoName, biomeName) {
  var eco = normalizar(ecoName);
  var bio = normalizar(biomeName);

  if (eco.indexOf('caatinga enclaves') >= 0) { return 'MATA_ATLANTICA'; }
  if (eco.indexOf('cerrado') >= 0) { return 'CERRADO'; }
  if (eco.indexOf('pantanal') >= 0) { return 'PANTANAL'; }
  if (eco.indexOf('uruguayan savanna') >= 0) { return 'PAMPA'; }
  if (eco.indexOf('caatinga') >= 0) { return 'CAATINGA'; }

  for (var i = 0; i < ECO_MATA_ATLANTICA.length; i++) {
    if (eco.indexOf(ECO_MATA_ATLANTICA[i]) >= 0) { return 'MATA_ATLANTICA'; }
  }

  // Floresta umida com nome desconhecido pode ser Amazonia ou Mata Atlantica:
  // deixa como nao conclusivo em vez de chutar.
  if (bio.indexOf('moist broadleaf') >= 0) { return null; }
  if (bio.indexOf('grasslands, savannas') >= 0 && bio.indexOf('temperate') >= 0) { return 'PAMPA'; }
  if (bio.indexOf('flooded grasslands') >= 0) { return 'PANTANAL'; }
  if (bio.indexOf('dry broadleaf') >= 0 || bio.indexOf('xeric') >= 0) { return 'CAATINGA'; }
  if (bio.indexOf('grasslands, savannas') >= 0) { return 'CERRADO'; }
  return null;
}

function detectarBioma(aoi, callback) {
  var eco = ee.FeatureCollection('RESOLVE/ECOREGIONS/2017').filterBounds(aoi);
  var comArea = eco.map(function(f) {
    var inter = ee.Feature(f).geometry().intersection(aoi, ee.ErrorMargin(100));
    return ee.Feature(f).set('area_int', inter.area(100));
  });
  var n = comArea.size();
  var maior = ee.Feature(ee.Algorithms.If(
    n.gt(0),
    comArea.sort('area_int', false).first(),
    ee.Feature(null, {ECO_NAME: null, BIOME_NAME: null})
  ));
  ee.Dictionary({
    eco: maior.get('ECO_NAME'),
    wwf: maior.get('BIOME_NAME'),
    n:   n
  }).evaluate(callback);
}

// SIRGAS 2000 / UTM do fuso do centroide (sul: 31960 + fuso; norte: 31954 + fuso).
function epsgSirgasUtm(lon, lat) {
  var fuso = Math.floor((lon + 180) / 6) + 1;
  return 'EPSG:' + String(lat < 0 ? 31960 + fuso : 31954 + fuso);
}


// ---------------------------------------------------------------------------
// Processamento
// ---------------------------------------------------------------------------
function mascarar(image) {
  var img = ee.Image(image);   // dentro do map as vezes chega como Feature
  return img.updateMask(img.select('cs_cdf').gte(CS_LIMIAR))
            .select(BANDAS)
            .divide(10000)
            .copyProperties(img, ['system:time_start']);
}

// filterDate exclui a data final, entao a janela vai ate o dia 1 do mes
// seguinte ao ultimo mes.
function janelaDoAno(P, ano) {
  var cruza = (P.mes_fim < P.mes_ini);
  var a = ee.Number(ano);
  var ini = ee.Date.fromYMD(a, P.mes_ini, 1);
  var fim = cruza
    ? ee.Date.fromYMD(a.add(1), P.mes_fim, 1).advance(1, 'month')
    : ee.Date.fromYMD(a, P.mes_fim, 1).advance(1, 'month');
  return {ini: ini, fim: fim};
}

// So entram anos com a janela completa.
function janelaFechada(P, ano) {
  var cruza = (P.mes_fim < P.mes_ini);
  var anoFim = cruza ? ano + 1 : ano;
  var fim = new Date(Date.UTC(anoFim, P.mes_fim, 1));  // mes em JS comeca do 0
  return new Date() >= fim;
}

function fracoes(image, P) {
  var img = ee.Image(image);
  var f = img.unmix([P.em_veg, P.em_solo, P.em_sombra], true, true)
             .rename(['f_veg','f_solo','f_sombra']);
  var ndvi = img.normalizedDifference(['B8','B4']).rename('ndvi');
  return img.addBands(f).addBands(ndvi);
}

function montarProdutos(aoi, P, anos) {
  var regiao = aoi.buffer(MARGEM);

  var anual = ee.ImageCollection(ee.List(anos).map(function(ano) {
    var j = janelaDoAno(P, ano);
    var col = ee.ImageCollection(S2.filterDate(j.ini, j.fim)
                .filterBounds(regiao)
                .linkCollection(CS, ['cs_cdf'])
                .map(mascarar));
    // n_obs = quantas imagens validas entraram na mediana de cada pixel
    var comp = fracoes(col.median().clip(regiao), P)
                 .addBands(col.select('B4').count().rename('n_obs').clip(regiao));
    return comp.set({ano: ano, n_cenas: col.size(), 'system:time_start': j.ini.millis()});
  }));

  var degAnual = anual.map(function(image) {
    var img = ee.Image(image);
    var d = img.select('f_solo').gte(P.lim_solo)
               .and(img.select('ndvi').lte(P.lim_ndvi))
               .rename('classificado');
    return d.set('ano', img.get('ano'))
            .set('system:time_start', img.get('system:time_start'));
  });

  var anoIni = anos[0], anoFim = anos[anos.length - 1];
  var imgIni = ee.Image(anual.filter(ee.Filter.eq('ano', anoIni)).first());
  var imgFim = ee.Image(anual.filter(ee.Filter.eq('ano', anoFim)).first());

  // dNdvi = final - inicial. A reducao e o mesmo valor com sinal trocado,
  // mantida so onde o NDVI caiu, sem classes.
  var dNdvi = imgFim.select('ndvi').subtract(imgIni.select('ndvi')).rename('d_ndvi');
  var reducaoNdvi = dNdvi.multiply(-1)
    .rename('reducao_ndvi')
    .updateMask(dNdvi.lt(0))
    .clip(aoi);

  return {
    aoi: aoi, regiao: regiao, anos: anos,
    anoIni: anoIni, anoFim: anoFim,
    anual: anual, degAnual: degAnual,
    persistencia: degAnual.sum().rename('n_anos').clip(aoi),
    reducaoNdvi: reducaoNdvi,
    dNdvi: dNdvi.clip(aoi),
    classificadoFim: ee.Image(degAnual.filter(ee.Filter.eq('ano', anoFim)).first()).clip(aoi)
  };
}

function rgbDoAno(geom, P, ano) {
  var j = janelaDoAno(P, ano);
  return ee.ImageCollection(
      S2.filterDate(j.ini, j.fim).filterBounds(geom)
        .linkCollection(CS, ['cs_cdf'])
        .map(function(image) {
          var img = ee.Image(image);
          return img.updateMask(img.select('cs_cdf').gte(CS_LIMIAR))
                    .select(['B4','B3','B2']).divide(10000);
        })
    ).median().clip(geom);
}

// Area (ha) dos pixels iguais a 1. Com crs e o numero de reporte, na mesma
// grade das exportacoes; sem crs e so o numero de tela.
function areaHa(mascaraBinaria, geom, escala, crs) {
  var args = {reducer: ee.Reducer.sum(), geometry: geom, scale: escala,
              maxPixels: 1e13};
  if (crs) { args.crs = crs; args.tileScale = 4; } else { args.bestEffort = true; }
  return ee.Number(mascaraBinaria.selfMask()
    .multiply(ee.Image.pixelArea()).divide(1e4)
    .reduceRegion(args).values().get(0));
}

function areaTotalHa(geom, escala, crs) {
  var args = {reducer: ee.Reducer.sum(), geometry: geom, scale: escala,
              maxPixels: 1e13};
  if (crs) { args.crs = crs; }
  return ee.Number(ee.Image.pixelArea().divide(1e4).reduceRegion(args).get('area'));
}

function areasPorClasse(img, geom, escala, crs, tela) {
  var args = {reducer: ee.Reducer.sum().group({groupField: 1, groupName: 'classe'}),
              geometry: geom, scale: escala, maxPixels: 1e13};
  if (crs) { args.crs = crs; }
  if (tela) { args.bestEffort = true; } else { args.tileScale = 4; }
  return ee.List(ee.Image.pixelArea().divide(1e4)
    .addBands(img.toInt().rename('classe'))
    .reduceRegion(args).get('groups'));
}

function mediaNaRegiao(img, geom, escala, crs) {
  var args = {reducer: ee.Reducer.mean(), geometry: geom, scale: escala,
              maxPixels: 1e13, tileScale: 4};
  if (crs) { args.crs = crs; }
  return img.reduceRegion(args);
}

// Chuva acumulada na janela (CHIRPS, pixel de ~5,5 km). Se a area for menor
// que um pixel, a media volta vazia e usamos o valor no centroide.
function chuvaNaJanela(P, ano, geom) {
  var j = janelaDoAno(P, ano);
  var soma = CHIRPS.filterDate(j.ini, j.fim).select('precipitation').sum();
  var media = soma.reduceRegion({reducer: ee.Reducer.mean(), geometry: geom,
                                 scale: 5566, maxPixels: 1e9, bestEffort: true})
                  .get('precipitation');
  var centro = soma.reduceRegion({reducer: ee.Reducer.first(),
                                  geometry: geom.centroid(100), scale: 5566})
                   .get('precipitation');
  return ee.Algorithms.If(media, media, centro);
}


// ---------------------------------------------------------------------------
// Estado e interface
// ---------------------------------------------------------------------------
var estado = {aoi: null, produtos: null, origemAoi: null, epsg: null,
              biomaVerificado: false, naMataAtlantica: null, ecorregiao: null,
              idExecucao: null};

var CINZA = '#5f6b63', VERDE = '#1b3a2b', VERMELHO = '#b3261e';

function titulo(txt) {
  return ui.Label(txt, {fontWeight: 'bold', fontSize: '14px', color: VERDE,
                        margin: '12px 0 4px 0'});
}
function subtitulo(txt) {
  return ui.Label(txt, {fontWeight: 'bold', fontSize: '12px', color: VERDE,
                        margin: '8px 0 2px 0'});
}
function nota(txt) {
  return ui.Label(txt, {fontSize: '11px', color: CINZA, margin: '2px 0 6px 0'});
}

var painel = ui.Panel({style: {width: '400px', padding: '10px'}});

painel.add(ui.Label('Mudanca de cobertura', {
  fontWeight: 'bold', fontSize: '18px', color: VERDE, margin: '0 0 2px 0'}));
painel.add(ui.Label('Propriedades militares  |  Mata Atlantica  |  Sentinel-2  |  ' +
                    ANO_INI + '-' + ANO_FIM, {
  fontSize: '12px', color: CINZA, margin: '0 0 8px 0'}));

// Passo 1: area
painel.add(titulo('1. Area de interesse'));
painel.add(nota('Desenhe um poligono no mapa da esquerda, ou informe um asset.'));

var txtAsset = ui.Textbox({placeholder: 'projects/.../assets/limite_om',
                           style: {width: '250px'}});
var btnAsset = ui.Button({label: 'Usar asset', style: {stretch: 'horizontal'}});
painel.add(ui.Panel([txtAsset, btnAsset], ui.Panel.Layout.flow('horizontal')));

var lblArea = ui.Label('Nenhuma area definida.', {fontSize: '12px', margin: '4px 0'});
painel.add(lblArea);

// Id da execucao (sigla_data_hora), usado nos nomes dos arquivos exportados.
function gerarIdExecucao(P) {
  var d = new Date();
  function z(n) { return (n < 10 ? '0' : '') + n; }
  return P.sigla + '_' +
    d.getFullYear() + z(d.getMonth() + 1) + z(d.getDate()) + '_' +
    z(d.getHours()) + z(d.getMinutes()) + z(d.getSeconds());
}

function idAreaAtual(P) {
  return estado.idExecucao || P.sigla;
}

// Passo 2: bioma
painel.add(titulo('2. Bioma'));
var btnBioma = ui.Button({label: 'Verificar bioma', style: {stretch: 'horizontal'}});
var lblEco = ui.Label('', {fontSize: '11px', color: CINZA, margin: '4px 0'});
var chkFora = ui.Checkbox({
  label: 'Processar mesmo assim, ciente de que o resultado nao tem verificacao',
  value: false, style: {shown: false, fontSize: '11px', color: VERMELHO}});
painel.add(btnBioma);
painel.add(lblEco);
painel.add(chkFora);
painel.add(nota('Os parametros foram testados apenas na Mata Atlantica.'));

// Passo 3: processar
painel.add(titulo('3. Processar'));
var btnProcessar = ui.Button({label: 'Processar serie ' + ANO_INI + '-' + ANO_FIM,
                              style: {stretch: 'horizontal'}});
painel.add(btnProcessar);
var lblStatus = ui.Label('', {fontSize: '12px', margin: '4px 0'});
painel.add(lblStatus);

var painelParam = ui.Panel({style: {shown: false, margin: '4px 0'}});
painel.add(painelParam);
var painelResultado = ui.Panel({style: {shown: false}});
painel.add(painelResultado);
var painelPixel = ui.Panel({style: {shown: false}});
painel.add(painelPixel);
var painelExport = ui.Panel({style: {shown: false}});
painel.add(painelExport);

// Mapas lado a lado com cortina: esquerda = ano inicial, direita = ano final.
var mapaA = ui.Map(), mapaB = ui.Map();
mapaA.setOptions('SATELLITE');
mapaB.setOptions('SATELLITE');
mapaA.setControlVisibility({drawingToolsControl: true});
mapaB.setControlVisibility({drawingToolsControl: false});
ui.Map.Linker([mapaA, mapaB]);

var split = ui.SplitPanel({firstPanel: mapaA, secondPanel: mapaB,
                           orientation: 'horizontal', wipe: true,
                           style: {stretch: 'both'}});

ui.root.clear();
ui.root.setLayout(ui.Panel.Layout.flow('horizontal'));
ui.root.add(painel);
ui.root.add(ui.Panel({widgets: [split], style: {stretch: 'both'}}));

var rotuloA = ui.Label(String(ANO_INI), {position: 'top-center', fontWeight: 'bold'});
var rotuloB = ui.Label(String(ANO_FIM), {position: 'top-center', fontWeight: 'bold'});
mapaA.add(rotuloA);
mapaB.add(rotuloB);

var desenho = mapaA.drawingTools();
desenho.setShape('polygon');
desenho.setDrawModes(['polygon', 'rectangle']);
desenho.setLinked(false);


// ---------------------------------------------------------------------------
// Legendas
// ---------------------------------------------------------------------------
function caixaLegenda(cor, texto) {
  return ui.Panel([
    ui.Label('', {backgroundColor: cor, padding: '8px', margin: '0 6px 2px 0',
                  border: '1px solid #999'}),
    ui.Label(texto, {fontSize: '11px', margin: '0 0 2px 0'})
  ], ui.Panel.Layout.flow('horizontal'));
}

var PAL_DNDVI = ['ffffcc', 'ffeda0', 'fed976', 'feb24c', 'fd8d3c',
                 'fc4e2a', 'e31a1c', 'bd0026', '800026'];
var PAL_PERS = ['ffffcc', 'fd8d3c', '800026'];

function estatisticasReducaoNdvi(prod) {
  var args = {reducer: ee.Reducer.minMax(), geometry: prod.aoi, scale: ESCALA,
              maxPixels: 1e13, tileScale: 4};
  if (estado.epsg) { args.crs = estado.epsg; }
  return prod.reducaoNdvi.reduceRegion(args);
}

// A rampa vai do minimo ao maximo da area; as cores nao sao classes.
function legendaReducaoNdvi(minimo, maximo) {
  var meio = (minimo + maximo) / 2;
  var p = ui.Panel({style: {stretch: 'horizontal', margin: '2px 0 6px 0'}});

  p.add(ui.Label('Reducao de NDVI (NDVI inicial - NDVI final)', {
    fontSize: '11px', margin: '0 0 2px 0'}));
  p.add(ui.Thumbnail({
    image: ee.Image.pixelLonLat().select('longitude'),
    params: {bbox: [0, 0, 1, 0.1], dimensions: '260x18', format: 'png',
             min: 0, max: 1, palette: PAL_DNDVI},
    style: {stretch: 'horizontal', margin: '0 0 2px 0', maxHeight: '18px'}
  }));
  p.add(ui.Panel([
    ui.Label(minimo.toFixed(3), {fontSize: '10px', width: '33%', textAlign: 'left'}),
    ui.Label(meio.toFixed(3),   {fontSize: '10px', width: '34%', textAlign: 'center'}),
    ui.Label(maximo.toFixed(3), {fontSize: '10px', width: '33%', textAlign: 'right'})
  ], ui.Panel.Layout.flow('horizontal'), {stretch: 'horizontal'}));
  p.add(nota('Escala ajustada a cada area. As cores nao indicam classes de severidade.'));
  return p;
}


// ---------------------------------------------------------------------------
// Acoes
// ---------------------------------------------------------------------------
function definirAoi(geom, origem) {
  estado.aoi = geom;
  estado.origemAoi = origem;
  estado.produtos = null;
  estado.epsg = null;
  estado.biomaVerificado = false;
  estado.naMataAtlantica = null;
  estado.ecorregiao = null;
  estado.idExecucao = null;
  lblEco.setValue('');
  chkFora.setValue(false);
  chkFora.style().set('shown', false);
  painelResultado.style().set('shown', false);
  painelExport.style().set('shown', false);
  painelPixel.style().set('shown', false);

  geom.area(100).divide(1e4).evaluate(function(ha, err) {
    if (err) { lblArea.setValue('Erro ao ler a area: ' + err); return; }
    lblArea.setValue('Area definida (' + origem + '): ' + ha.toFixed(2) + ' ha');
  });
  // Map.centerObject as vezes falha com erro de maxError; setCenter nao.
  // O centroide tambem define o fuso UTM das exportacoes.
  geom.centroid(10).coordinates().evaluate(function(c) {
    if (c) {
      mapaA.setCenter(c[0], c[1], 12);
      estado.epsg = epsgSirgasUtm(c[0], c[1]);
    }
  });
}

desenho.onDraw(function(geom) {
  // fica so o ultimo poligono desenhado
  while (desenho.layers().length() > 1) { desenho.layers().remove(desenho.layers().get(0)); }
  definirAoi(ee.Geometry(geom), 'desenho');
});
desenho.onEdit(function(geom) { definirAoi(ee.Geometry(geom), 'desenho'); });

btnAsset.onClick(function() {
  var id = txtAsset.getValue();
  if (!id) { lblArea.setValue('Informe o id do asset.'); return; }
  lblArea.setValue('Carregando asset...');
  var g = ee.FeatureCollection(id).geometry();
  g.area(100).evaluate(function(a, err) {
    if (err || a === null) {
      lblArea.setValue('Nao foi possivel abrir o asset. Confira o id e a permissao.');
      return;
    }
    definirAoi(g, 'asset');
  });
});

btnBioma.onClick(function() {
  if (!estado.aoi) { lblEco.setValue('Defina a area primeiro.'); return; }
  lblEco.setValue('Verificando...');
  detectarBioma(estado.aoi, function(res, err) {
    estado.biomaVerificado = true;
    if (err || !res || res.n === 0 || !res.eco) {
      estado.naMataAtlantica = null;
      lblEco.setValue('Nao foi possivel identificar a ecorregiao. Confira no mapa ' +
                      'do IBGE se a area esta na Mata Atlantica.');
      lblEco.style().set('color', VERMELHO);
      chkFora.style().set('shown', true);
      return;
    }
    var chave = biomaDaEcorregiao(res.eco, res.wwf);
    estado.ecorregiao = res.eco;
    estado.naMataAtlantica = (chave === 'MATA_ATLANTICA') ? true : (chave ? false : null);
    if (estado.naMataAtlantica === true) {
      lblEco.setValue('Mata Atlantica (ecorregiao predominante: ' + res.eco + ').');
      lblEco.style().set('color', CINZA);
      chkFora.setValue(false);
      chkFora.style().set('shown', false);
    } else if (estado.naMataAtlantica === false) {
      lblEco.setValue('Area FORA da Mata Atlantica: ecorregiao predominante ' + res.eco +
                      ' (' + ROTULO_BIOMA[chave] + '). Os parametros so foram testados ' +
                      'na Mata Atlantica; aqui o resultado nao tem verificacao.');
      lblEco.style().set('color', VERMELHO);
      chkFora.style().set('shown', true);
    } else {
      lblEco.setValue('Ecorregiao predominante: ' + res.eco + '. Nao foi possivel ' +
                      'confirmar o bioma pelo nome; confira no mapa do IBGE.');
      lblEco.style().set('color', VERMELHO);
      chkFora.style().set('shown', true);
    }
  });
});

function mostrarParametros(P, anos) {
  painelParam.clear();
  painelParam.style().set('shown', true);
  painelParam.add(titulo('Parametros usados'));
  if (estado.naMataAtlantica !== true) {
    painelParam.add(ui.Label('ATENCAO: area fora da Mata Atlantica ou nao confirmada. ' +
                             'Resultado sem verificacao.',
                             {fontSize: '11px', color: VERMELHO, fontWeight: 'bold'}));
  }
  painelParam.add(ui.Label(
    'Bioma: ' + P.rotulo + '  |  janela: mes ' + P.mes_ini + ' a ' + P.mes_fim + '\n' +
    'Limiares: f_solo >= ' + P.lim_solo + ' e NDVI <= ' + P.lim_ndvi + '\n' +
    'Cloud Score+ >= ' + CS_LIMIAR + '  |  anos: ' + anos.join(', ') + '\n' +
    'Exportacao: ' + estado.epsg + ', ' + ESCALA + ' m\n' +
    'Endmembers (pixels puros da propria imagem):\n' +
    '  veg ' + P.em_veg.join(', ') + '\n' +
    '  solo ' + P.em_solo.join(', ') + '\n' +
    '  sombra ' + P.em_sombra.join(', ') + '\n' +
    P.nota,
    {fontSize: '11px', color: CINZA, margin: '2px 0 6px 0', whiteSpace: 'pre'}));
}

btnProcessar.onClick(function() {
  if (!estado.aoi) { lblStatus.setValue('Defina a area primeiro.'); return; }
  if (!estado.epsg) {
    lblStatus.setValue('A area ainda esta sendo lida. Aguarde alguns segundos e tente de novo.');
    return;
  }
  if (!estado.biomaVerificado) { lblStatus.setValue('Verifique o bioma (passo 2).'); return; }
  if (estado.naMataAtlantica !== true && !chkFora.getValue()) {
    lblStatus.setValue('Area fora da Mata Atlantica ou nao confirmada: marque a opcao ' +
                       'do passo 2 para processar sem verificacao.');
    return;
  }

  var P = PARAM;
  estado.idExecucao = gerarIdExecucao(P);

  var anos = [], fora = [];
  for (var a = ANO_INI; a <= ANO_FIM; a++) {
    if (janelaFechada(P, a)) { anos.push(a); } else { fora.push(a); }
  }
  if (anos.length < 2) {
    lblStatus.setValue('Ainda nao ha anos suficientes com a janela completa.');
    return;
  }

  lblStatus.setValue('Processando ' + anos.length + ' anos...' +
    (fora.length ? ' (ficaram de fora, janela incompleta: ' + fora.join(', ') + ')' : ''));

  var prod = montarProdutos(estado.aoi, P, anos);
  estado.produtos = prod;
  mostrarParametros(P, anos);
  desenharMapas(prod, P);
  calcularResultados(prod, P);
  montarSaidas(prod, P);

  // A amostra vai para a aba Tasks a cada processamento; rode so a da
  // execucao que for usada no trabalho.
  exportarAmostra(prod, P, idAreaAtual(P), N_VALIDACAO_POR_ESTRATO);
  lblStatus.setValue(
    'Pronto. Amostra de validacao (' + N_VALIDACAO_POR_ESTRATO +
    ' pontos por estrato) e areas dos estratos enviadas para a aba Tasks.');
});

function desenharMapas(prod, P) {
  mapaA.layers().reset();
  mapaB.layers().reset();

  var visRGB = {min: 0.02, max: 0.30};
  mapaA.addLayer(rgbDoAno(prod.aoi, P, prod.anoIni), visRGB, 'RGB ' + prod.anoIni);
  mapaB.addLayer(rgbDoAno(prod.aoi, P, prod.anoFim), visRGB, 'RGB ' + prod.anoFim);

  var clsIni = ee.Image(prod.degAnual.filter(ee.Filter.eq('ano', prod.anoIni)).first()).clip(prod.aoi);
  mapaA.addLayer(clsIni.selfMask(), {palette: ['e34a33']}, 'Classificado ' + prod.anoIni);
  mapaB.addLayer(prod.classificadoFim.selfMask(), {palette: ['e34a33']}, 'Classificado ' + prod.anoFim);

  mapaB.addLayer(prod.persistencia.selfMask(),
    {min: 1, max: prod.anos.length, palette: PAL_PERS}, 'Persistencia', false);

  // A camada entra ja na ordem certa; a rampa e ajustada quando o min/max chegar.
  var camadaReducao = ui.Map.Layer(prod.reducaoNdvi, {min: 0, max: 1, palette: PAL_DNDVI},
                                   'Reducao de NDVI', false);
  mapaB.layers().add(camadaReducao);

  estatisticasReducaoNdvi(prod).evaluate(function(stats, err) {
    if (err || !stats || stats.reducao_ndvi_min === null ||
        stats.reducao_ndvi_max === null) {
      camadaReducao.setName('Reducao de NDVI (sem pixels com reducao)');
      return;
    }
    var minimo = stats.reducao_ndvi_min;
    var maximo = stats.reducao_ndvi_max;
    var maxVis = (maximo > minimo) ? maximo : minimo + 1e-6;
    camadaReducao.setVisParams({min: minimo, max: maxVis, palette: PAL_DNDVI});
    camadaReducao.setName('Reducao de NDVI [' + minimo.toFixed(3) + ' a ' +
                          maximo.toFixed(3) + ']');
  });

  var contorno = ee.Image().paint(ee.FeatureCollection([ee.Feature(prod.aoi)]), 0, 2);
  mapaA.addLayer(contorno, {palette: '00ffff'}, 'Limite da propriedade');
  mapaB.addLayer(contorno, {palette: '00ffff'}, 'Limite da propriedade');

  rotuloA.setValue(String(prod.anoIni));
  rotuloB.setValue(String(prod.anoFim));
}

function calcularResultados(prod, P) {
  painelResultado.clear();
  painelResultado.style().set('shown', true);
  painelResultado.add(titulo('Resultados dentro do poligono'));
  var lbl = ui.Label('Calculando...', {fontSize: '12px'});
  painelResultado.add(lbl);

  var totalHa = ee.Number(ee.Image.pixelArea().divide(1e4).reduceRegion({
    reducer: ee.Reducer.sum(), geometry: prod.aoi,
    scale: ESCALA_UI, maxPixels: 1e13, bestEffort: true}).values().get(0));

  var clsHa = areaHa(prod.classificadoFim, prod.aoi, ESCALA_UI);
  var redHa = areaHa(prod.dNdvi.lt(0), prod.aoi, ESCALA_UI);
  var statsNdvi = estatisticasReducaoNdvi(prod);

  var painelLegendaNdvi = ui.Panel();

  ee.Dictionary({
    total: totalHa,
    cls: clsHa,
    red: redHa,
    ndviMin: statsNdvi.get('reducao_ndvi_min'),
    ndviMax: statsNdvi.get('reducao_ndvi_max')
  }).evaluate(function(r, err) {
    if (err) { lbl.setValue('Erro no calculo: ' + err); return; }
    var pct = r.total ? (100 * r.cls / r.total) : 0;

    lbl.setValue(
      'Area total: ' + r.total.toFixed(2) + ' ha\n' +
      'Classificada em ' + prod.anoFim + ': ' + r.cls.toFixed(2) + ' ha (' +
        pct.toFixed(2) + '%)\n' +
      'Area com reducao de NDVI: ' + r.red.toFixed(2) + ' ha' +
      ((r.ndviMin !== null && r.ndviMax !== null)
        ? '\nReducao de NDVI: min ' + r.ndviMin.toFixed(3) +
          ' | max ' + r.ndviMax.toFixed(3)
        : '\nNenhum pixel com reducao de NDVI.'));
    lbl.style().set('whiteSpace', 'pre');

    painelLegendaNdvi.clear();
    if (r.ndviMin !== null && r.ndviMax !== null) {
      painelLegendaNdvi.add(legendaReducaoNdvi(r.ndviMin, r.ndviMax));
    } else {
      painelLegendaNdvi.add(nota('Reducao de NDVI: nenhum pixel com queda do indice.'));
    }
  });

  painelResultado.add(nota('Numeros de tela, calculados a ' + ESCALA_UI + ' m. ' +
    'Para o relatorio, use as tabelas exportadas (' + ESCALA + ' m).'));

  painelResultado.add(titulo('Legenda'));
  painelResultado.add(ui.Label('Persistencia (anos classificado)', {fontSize: '11px'}));
  painelResultado.add(caixaLegenda(PAL_PERS[0], '1'));
  painelResultado.add(caixaLegenda(PAL_PERS[1], 'intermediaria'));
  painelResultado.add(caixaLegenda(PAL_PERS[2], String(prod.anos.length) + ' (todos os anos)'));
  painelResultado.add(painelLegendaNdvi);

  painelResultado.add(nota('Clique em um ponto do mapa para ver a serie daquele pixel.'));
}


// ---------------------------------------------------------------------------
// Serie do pixel clicado
// ---------------------------------------------------------------------------
function serieDoPonto(prod, P, ponto) {
  return ee.FeatureCollection(prod.anual.map(function(image) {
    var img = ee.Image(image);
    var v = img.select(['f_solo','ndvi']).reduceRegion({
      reducer: ee.Reducer.first(), geometry: ponto, scale: ESCALA});
    return ee.Feature(null, {ano: img.get('ano'),
                             f_solo: v.get('f_solo'),
                             ndvi: v.get('ndvi')});
  }));
}

function aoClicar(coords) {
  if (!estado.produtos) { return; }
  var prod = estado.produtos, P = PARAM;
  var ponto = ee.Geometry.Point([coords.lon, coords.lat]);

  painelPixel.clear();
  painelPixel.style().set('shown', true);
  painelPixel.add(titulo('Serie do pixel'));
  var lblP = ui.Label('Lendo...', {fontSize: '12px'});
  painelPixel.add(lblP);

  serieDoPonto(prod, P, ponto).evaluate(function(fc, err) {
    painelPixel.clear();
    painelPixel.add(titulo('Serie do pixel'));
    if (err || !fc || !fc.features.length) {
      painelPixel.add(ui.Label('Sem dado neste ponto.', {fontSize: '12px'}));
      return;
    }

    var tabela = [['Ano', 'Fracao de solo', 'NDVI']];
    var classificados = [];
    for (var i = 0; i < fc.features.length; i++) {
      var p = fc.features[i].properties;
      if (p.f_solo === null || p.ndvi === null) { continue; }
      tabela.push([String(p.ano), p.f_solo, p.ndvi]);
      if (p.f_solo >= P.lim_solo && p.ndvi <= P.lim_ndvi) {
        classificados.push(p.ano);
      }
    }

    if (tabela.length < 2) {
      painelPixel.add(ui.Label('Sem observacao valida neste ponto.', {fontSize: '12px'}));
      return;
    }

    painelPixel.add(ui.Chart(tabela, 'LineChart', {
      title: 'Serie temporal do pixel',
      hAxis: {title: 'Ano'},
      vAxis: {viewWindow: {min: 0, max: 1}},
      series: {0: {color: 'b35806'}, 1: {color: '1b3a2b'}},
      legend: {position: 'bottom'},
      height: 220
    }));
    painelPixel.add(ui.Label(
      classificados.length
        ? 'Classificado em: ' + classificados.join(', ') +
          '  (' + classificados.length + ' de ' + prod.anos.length + ' anos)'
        : 'Nunca classificado na serie.',
      {fontSize: '12px', margin: '2px 0'}));
  });
}

mapaA.onClick(aoClicar);
mapaB.onClick(aoClicar);


// ---------------------------------------------------------------------------
// Amostra de validacao e exportacoes
// ---------------------------------------------------------------------------
// Amostra estratificada pelo mapa do ultimo ano (classe_alg 1 = classificado,
// 0 = nao classificado). O SHP leva a classe do mapa: para interpretar sem
// ver essa classe, gere o KML so com o id antes de abrir no Google Earth Pro.
function exportarAmostra(prod, P, idArea, nPorEstrato) {
  var img = prod.classificadoFim.rename('classe_alg');

  var args = {numPoints: nPorEstrato, classBand: 'classe_alg', region: prod.aoi,
              scale: ESCALA, seed: SEMENTE_AMOSTRA, geometries: true, tileScale: 4};
  if (estado.epsg) { args.projection = ee.Projection(estado.epsg); }

  var pontos = img.stratifiedSample(args)
                  .randomColumn('ordem', SEMENTE_ORDEM).sort('ordem');

  var lista = pontos.toList(pontos.size());
  var comId = ee.FeatureCollection(
    ee.List.sequence(0, pontos.size().subtract(1)).map(function(j) {
      var f = ee.Feature(lista.get(j));
      return ee.Feature(f.geometry(), {
        id: ee.String(idArea + '-').cat(ee.Number(j).add(1).int().format('%03d')),
        classe_alg: f.get('classe_alg')
      });
    }));

  Export.table.toDrive({
    collection: comId,
    description: 'amostras_' + idArea,
    folder: PASTA_EXPORT,
    fileFormat: 'SHP',
    selectors: ['id', 'classe_alg']
  });

  // Areas dos estratos, para os pesos do estimador de acuracia.
  var crs = estado.epsg;
  var estratos = ee.FeatureCollection([ee.Feature(null, {
    identificador: idArea, bioma: P.rotulo, ano_mapa: prod.anoFim,
    area_estrato1_ha: areaHa(prod.classificadoFim.eq(1), prod.aoi, ESCALA, crs),
    area_estrato0_ha: areaHa(prod.classificadoFim.eq(0), prod.aoi, ESCALA, crs),
    pontos_por_estrato: nPorEstrato, semente: SEMENTE_AMOSTRA,
    epsg: crs, versao_app: VERSAO
  })]);
  Export.table.toDrive({
    collection: estratos, description: 'estratos_' + idArea,
    folder: PASTA_EXPORT, fileFormat: 'CSV',
    selectors: ['identificador', 'bioma', 'ano_mapa', 'area_estrato1_ha',
                'area_estrato0_ha', 'pontos_por_estrato', 'semente', 'epsg',
                'versao_app']});
}

function montarSaidas(prod, P) {
  painelExport.clear();
  painelExport.style().set('shown', true);
  painelExport.add(titulo('4. Saidas'));
  painelExport.add(nota('Cada botao manda uma tarefa para a aba Tasks; clique em ' +
    'RUN la para gravar no Drive (pasta ' + PASTA_EXPORT + '). Os arquivos levam o ' +
    'id da execucao no nome.'));

  function idArea() { return idAreaAtual(P); }

  painelExport.add(subtitulo('Exportacoes'));

  function expImagem(img, nome, regiao) {
    Export.image.toDrive({image: img, description: nome, folder: PASTA_EXPORT,
                          region: regiao || prod.aoi, scale: ESCALA,
                          crs: estado.epsg, maxPixels: 1e13});
  }

  // Rasters de classe: fora da area (e onde nao houve imagem valida) sai 255,
  // para nao se confundir com a classe 0.
  function expClasse(img, nome) {
    expImagem(img.unmask(SEM_DADO).toByte(), nome, prod.aoi.bounds());
  }

  function botao(rotulo, acao) {
    var b = ui.Button({label: rotulo, style: {stretch: 'horizontal'}});
    b.onClick(function() { acao(); b.setLabel(rotulo + '  (na fila)'); });
    painelExport.add(b);
  }

  botao('Area de estudo (SHP)', function() {
    var areaVet = ee.FeatureCollection([
      ee.Feature(prod.aoi, {
        id_area: idArea(),
        area_ha: prod.aoi.area(1).divide(1e4),
        bioma: P.rotulo,
        origem: estado.origemAoi || 'nao_def',
        epsg: estado.epsg,
        versao: VERSAO
      })
    ]);
    // nomes curtos por causa do limite de 10 caracteres do DBF
    Export.table.toDrive({
      collection: areaVet,
      description: 'area_estudo_' + idArea(),
      folder: PASTA_EXPORT,
      fileFormat: 'SHP',
      selectors: ['id_area', 'area_ha', 'bioma', 'origem', 'epsg', 'versao']
    });
  });

  botao('Classificado ' + prod.anoFim, function() {
    expClasse(prod.classificadoFim, 'classificado_' + idArea() + '_' + prod.anoFim);
  });
  botao('Persistencia', function() {
    expClasse(prod.persistencia, 'persistencia_' + idArea());
  });
  botao('Reducao de NDVI (continua)', function() {
    expImagem(prod.reducaoNdvi.toFloat(), 'reducao_ndvi_' + idArea());
  });

  // As seis bandas que entram no MLME, em reflectancia (0-1), um arquivo por
  // ano. B11 e B12 (20 m) saem reamostradas para a grade de 10 m.
  botao('Composicoes multiespectrais (B2,B3,B4,B8,B11,B12)', function() {
    for (var i = 0; i < prod.anos.length; i++) {
      var ano = prod.anos[i];
      var comp6 = ee.Image(prod.anual.filter(ee.Filter.eq('ano', ano)).first())
        .select(BANDAS)
        .toFloat();
      expImagem(comp6, 'multiespectral_6bandas_' + idArea() + '_' + ano, prod.regiao);
    }
  });

  botao('RGB ' + prod.anoIni + ' e ' + prod.anoFim + ' (para figuras)', function() {
    var anosRgb = [prod.anoIni, prod.anoFim];
    for (var i = 0; i < anosRgb.length; i++) {
      expImagem(rgbDoAno(prod.regiao, P, anosRgb[i]).multiply(10000).toUint16(),
                'rgb_' + idArea() + '_' + anosRgb[i], prod.regiao);
    }
  });

  // Uma linha por ano: area classificada e o que pode explicar a variacao
  // (cenas, observacoes validas, chuva, medias de f_solo e NDVI).
  botao('Serie anual (CSV)', function() {
    var crs = estado.epsg;
    var total = areaTotalHa(prod.aoi, ESCALA, crs);
    var serie = ee.FeatureCollection(prod.anos.map(function(ano) {
      var comp = ee.Image(prod.anual.filter(ee.Filter.eq('ano', ano)).first());
      var cls  = ee.Image(prod.degAnual.filter(ee.Filter.eq('ano', ano)).first());
      var medias = mediaNaRegiao(comp.select(['f_solo', 'ndvi']), prod.aoi, ESCALA, crs);
      var obs = comp.select('n_obs').reduceRegion({
        reducer: ee.Reducer.median(), geometry: prod.aoi, scale: ESCALA,
        crs: crs, maxPixels: 1e13, tileScale: 4}).get('n_obs');
      var ha = areaHa(cls, prod.aoi, ESCALA, crs);
      return ee.Feature(null, {
        ano: ano,
        n_cenas: comp.get('n_cenas'),
        obs_validas_mediana: obs,
        chuva_mm: chuvaNaJanela(P, ano, prod.aoi),
        area_ha_propriedade: ha,
        area_total_ha_propriedade: total,
        pct_propriedade: ha.divide(total).multiply(100),
        media_fsolo_propriedade: medias.get('f_solo'),
        media_ndvi_propriedade: medias.get('ndvi')
      });
    }));
    Export.table.toDrive({
      collection: serie, description: 'serie_anual_' + idArea(),
      folder: PASTA_EXPORT, fileFormat: 'CSV',
      selectors: ['ano', 'n_cenas', 'obs_validas_mediana', 'chuva_mm',
                  'area_ha_propriedade', 'area_total_ha_propriedade', 'pct_propriedade',
                  'media_fsolo_propriedade', 'media_ndvi_propriedade']});
  });

  // Area de cada classe do mapa final e da persistencia (tabelas do relatorio).
  botao('Areas por classe (CSV)', function() {
    var crs = estado.epsg;
    function linhas(img, nome) {
      return ee.FeatureCollection(
        areasPorClasse(img, prod.aoi, ESCALA, crs, false).map(function(g) {
          var d = ee.Dictionary(g);
          return ee.Feature(null, {produto: nome, classe: d.get('classe'),
                                   area_ha: d.get('sum')});
        }));
    }
    var tabela = linhas(prod.classificadoFim, 'classificado_' + prod.anoFim)
      .merge(linhas(prod.persistencia, 'persistencia_anos'));
    Export.table.toDrive({
      collection: tabela, description: 'areas_por_classe_' + idArea(),
      folder: PASTA_EXPORT, fileFormat: 'CSV',
      selectors: ['produto', 'classe', 'area_ha']});
  });

  // Registro do que foi usado em cada execucao.
  botao('Parametros usados (CSV)', function() {
    var params = ee.FeatureCollection([ee.Feature(null, {
      identificador: idArea(),
      bioma: P.rotulo,
      ano_inicial: prod.anoIni, ano_final: prod.anoFim,
      anos_processados: prod.anos.join(' '),
      mes_ini: P.mes_ini, mes_fim: P.mes_fim,
      limiar_solo: P.lim_solo, limiar_ndvi: P.lim_ndvi,
      cloud_score_min: CS_LIMIAR, escala_m: ESCALA,
      epsg: estado.epsg,
      em_veg: P.em_veg.join(' '), em_solo: P.em_solo.join(' '),
      em_sombra: P.em_sombra.join(' '),
      endmembers_origem: 'pixels puros da imagem',
      origem_area: estado.origemAoi,
      ecorregiao: estado.ecorregiao || 'nao identificada',
      area_na_mata_atlantica: estado.naMataAtlantica === true ? 'sim' :
        (estado.naMataAtlantica === false ? 'nao (sem verificacao)' : 'nao confirmado'),
      versao_app: VERSAO
    })]);
    Export.table.toDrive({collection: params, description: 'parametros_' + idArea(),
                          folder: PASTA_EXPORT, fileFormat: 'CSV'});
  });
}


// ---------------------------------------------------------------------------
// Inicio
// ---------------------------------------------------------------------------
mapaA.setCenter(-43.5, -22.85, 11);
painel.add(nota('Serie fixa em ' + ANO_INI + '-' + ANO_FIM + '. Anos com a janela ' +
                'jun-set ainda incompleta ficam de fora.'));
