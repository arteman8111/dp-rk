import * as param from "./const.js"
import rungekutta from "./rungekutta.js";

// НУ
const app = document.querySelector("#app")
const table = document.createElement("table");
table.classList.add('app__table')
const tr_heade = document.createElement('tr')
tr_heade.classList.add('tr__header')
const tr_head = `
<tr class="tr__header">
<td>t, с</td>
<td>m, кг</td>
<td>Vx, м/с</td>
<td>Vy, м/с</td>
<td>x, м</td>
<td>y, м </td>
<td style="padding-left: 20px">V, м/с </td>
<td style="padding-left: 20px">r, м</td>
<td style="padding-left: 40px">thet, град</td>
<td style="padding-left: 20px">THET, град</td>
<td style="padding-left: 20px">alfa, град</td>
<td>fi, град</td>
</tr>
`
tr_heade.insertAdjacentHTML('beforeend', tr_head)
app.append(table) 
table.insertAdjacentElement('beforeend', tr_heade)
let dt = param.step;
let t = param.t0;
let t_prev = t;
let el = [param.m0, param.vx0, param.vy0, param.x0, param.y0, param.v0, param.r0,0,0,0,0];
let el_prev = el.slice();
function tr_iter(t,el){
    return `
        <td class="app__ceil">${t.toFixed(5)}</td>
        <td class="app__ceil">${el[0].toFixed(4)}</td>
        <td class="app__ceil">${el[1].toFixed(4)}</td>
        <td class="app__ceil">${el[2].toFixed(4)}</td>
        <td class="app__ceil">${el[3].toFixed(4)}</td>
        <td class="app__ceil">${el[4].toFixed(4)}</td>
        <td class="app__ceil">${el[5].toFixed(4)}</td>
        <td class="app__ceil">${el[6].toFixed(4)}</td>
        <td class="app__ceil">${el[7].toFixed(4)}</td>
        <td class="app__ceil">${el[8].toFixed(4)}</td>
        <td class="app__ceil">${el[9].toFixed(4)}</td>
        <td class="app__ceil">${el[10].toFixed(4)}</td>
    `
}
function html(t, el){
    const tr = document.createElement('tr')
    tr.insertAdjacentHTML('beforeend', tr_iter(t, el))
    table.insertAdjacentElement('beforeend', tr)
}

function init() {
    html(t, el)
    while (
        el[5] < param.vk
        // t < 20
        ) {
        if (t + dt > param.t1 && t < param.t1) {
            dt = param.t1 - t;
            t += dt;
            rungekutta(el, dt, t);
            html(t, el)
            dt = param.step - dt;
            t += dt;
            rungekutta(el, dt, t);
            html(t, el);
            dt = param.step
        }
        if (t + dt > param.t2 && t < param.t2) {
            dt = param.t2 - t;
            t += dt;
            rungekutta(el, dt, t);
            html(t, el)
            dt = param.step - dt;
            t += dt;
            rungekutta(el, dt, t);
            html(t, el)
            dt = param.step
        }
        el_prev = el.slice();
        t_prev = t;
        rungekutta(el, dt, t);
        t += dt;
        if (el[5] < param.vk) {
            html(t, el);
        }
        if (el[5] - param.eps_v > param.vk) {
            el = el_prev.slice();
            t = t_prev;
            dt = dt / 10;
            t += dt;
            rungekutta(el, dt, t);
        }
    }
}
init()